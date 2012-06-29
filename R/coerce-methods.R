dataFrameFromRange <- function(range, object, frame=0L, range.index=1L){
	## to do: change to S4 method and do dispatch on class of range
	if(missing(frame)) frame <- 200e3
	if(is(range, "RangedDataCNV")){
		rm <- IRanges::findOverlaps(range, featureData(object), maxgap=frame) ## RangesMatching
		sample.index <- match(sampleNames(range), sampleNames(object))
	} else {
		frange <- oligoClasses::makeFeatureGRanges(object)
		rm <- IRanges::findOverlaps(range, frange, maxgap=frame)
		sample.index <- match(sampleNames(range), sampleNames(object))
	}
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
	##df$range <- rep(i, nrow(df))##mm.df$query
	##dfList[[i]] <- df
	##df$range <- range.index
	if(is(range, "RangedDataCNV")){
		df$range <- paste("[", range.index, "] chr", chromosome(range), ", ID: ", sampleNames(obj), sep="")
	}
	if(is(range, "GRanges")){
		df$range <- paste("[", range.index, "] ", chromosome(range), ", ID: ", sampleNames(obj), sep="")
	}
	return(df)
}
