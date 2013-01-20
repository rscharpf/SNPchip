setMethod("dataFrame", signature(range="GRanges", data="gSet"),
	  function(range, data, ...){
		  dataFrameFromRange(range=range, object=data, ...)
	  })

setMethod("dataFrame", signature(range="GRanges", data="SummarizedExperiment"),
	  function(range, data, ...){
		  dataFrameSummarizedExperiment(range=range, object=data, ...)
	  })



setMethod("coerce", signature(from="BafLrrSetList", to="SummarizedExperiment"),
	  function(from, to){
		  ##nms <- varLabels(from@featureDataList[[1]])
		  chrom <- rep(paste("chr", chromosome(from), sep=""),
			       elementLengths(from))
		  pos <- unlist(position(from))
		  is.snp <- unlist(lapply(featureDataList(from), isSnp))
		  ## stack the featureDataList to make featureData
		  ## make granges object from featureData
		  sl <- getSequenceLengths(genomeBuild(from))
		  sl <- sl[unique(chrom)]

		  seqinfo <- Seqinfo(seqnames=unique(chrom),
				     genome="hg18")
		  gr <- GRanges(chrom, IRanges(pos,pos), isSnp=is.snp,
				seqlengths=sl,
				seqinfo=seqinfo)
		  names(gr) <- unlist(featureNames(from))
		  rlist <- lrr(from)
		  blist <- baf(from)
		  isff <- is(rlist[[1]], "ff")
		  if(isff) require("ff")
		  ##if(is(rlist[[1]], "ff")
		  rl <- lapply(rlist, "[", drop=FALSE) ##function(x) x[, ,drop=FALSE])
		  bl <- lapply(blist, "[", drop=FALSE) ##function(x) x[, ,drop=FALSE])
		  r <- do.call("rbind", rl)
		  b <- do.call("rbind", bl)
		  ##rownames(r) <- rownames(b) <- unlist(featureNames(from))
		  colData <- DataFrame(pData(from))
		  rownames(colData) <- sampleNames(from)
		  se <- SummarizedExperiment(assays=SimpleList(lrr=r, baf=b),
					     rowData=gr,
					     colData=colData)
		  return(se)
	  })

dataFrameSummarizedExperiment <- function(range, object, ...){
	range <- range[sampleNames(range) %in% colnames(object), ]
	grl <- split(range, sampleNames(range))
	if("maxgap" %in% names(list(...))){
		min.gapwidth <- list(...)[["maxgap"]]
		grl2 <- reduce(grl, min.gapwidth=min.gapwidth)
	} else grl2 <- reduce(grl)
	col.index <- match(names(grl2), colnames(object))
	selist <- foreach(gr=grl2, j=col.index) %do% subsetByOverlaps(object[, j], gr, ...)
	x <- unlist(lapply(selist, start))
	r <- unlist(lapply(selist, lrr))/100
	b <- unlist(lapply(selist, baf))/1000
	is.snp <- unlist(lapply(selist, isSnp))
	gr <- unlist(grl2)
	interval <- rep(seq_along(gr), elementLengths(selist))
	chrom <- rep(chromosome(gr), elementLengths(selist))
	id <- rep(names(gr), elementLengths(selist))
	## an interval may contain multiple CNVs.
	interval <- paste(chromosome(gr), " interval ", interval, ", ID: ", id, sep="")
	df <- data.frame(x=x, lrr=r, baf=b,
			 id=id,
			 is.snp=is.snp,
			 interval=interval)
	return(df)
}

dataFrameFromRange <- function(range, object, frame=0L, range.index=1L){
	## to do: change to S4 method and do dispatch on class of range
	if(missing(frame)) frame <- 200e3
	if(is(range, "RangedDataCNV")){
		rm <- IRanges::findOverlaps(range, featureData(object), maxgap=frame) ## RangesMatching
	} else {
		frange <- oligoClasses::makeFeatureGRanges(object)
		rm <- IRanges::findOverlaps(range, frange, maxgap=frame)
	}
	if(length(sampleNames(range))==0) {
		sample.index <- seq_len(ncol(object))
	} else  sample.index <- match(sampleNames(range), sampleNames(object))
	if(any(is.na(sample.index))) stop("sampleNames in RangedData do not match sampleNames in ", class(data), " object")
	sample.index <- unique(sample.index)
	mm <- IRanges::as.matrix(rm)
	mm.df <- data.frame(mm)
	mm.df$featureNames <- Biobase::featureNames(object)[mm.df$subject]
	marker.index <- mm.df$subject
	obj <- object[marker.index, sample.index]
	mm.df$subject <- match(mm.df$featureNames, featureNames(obj))
	##
	## coersion to data.frame
	##
	df <- as(obj, "data.frame")
	if(!missing(range.index)){
		df$range <- paste("[", range.index, "] ", chromosome(range), ", ID: ", sampleNames(obj), sep="")
	} else df$range <- paste(chromosome(range), ", ID: ", sampleNames(obj), sep="")
	return(df)
}


