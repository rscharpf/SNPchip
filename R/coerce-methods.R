dataFrameFromRange <- function(range, object, frame, range.index){
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
	mm <- as.matrix(rm)
	mm.df <- data.frame(mm)
	mm.df$featureNames <- featureNames(object)[mm.df$subject]
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
	df$range <- paste("[", range.index, "] chr ", chromosome(range), ", ID: ", sampleNames(obj), sep="")
	return(df)
}
