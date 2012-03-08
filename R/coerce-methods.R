dataFrameFromRange <- function(range, object, frame, range.index){
	rm <- findOverlaps(range, featureData(object), maxgap=frame) ## RangesMatching
	mm <- as.matrix(rm)
	mm.df <- data.frame(mm)
	mm.df$featureNames <- featureNames(object)[mm.df$subject]
	marker.index <- mm.df$subject
	sample.index <- match(sampleNames(range), sampleNames(object))
	if(any(is.na(sample.index))) stop("sampleNames in RangedData do not match sampleNames in ", class(data), " object")
	sample.index <- unique(sample.index)
	obj <- object[marker.index, sample.index]
	mm.df$subject <- match(mm.df$featureNames, featureNames(obj))
	##
	## coersion to data.frame
	##
	df <- as(obj, "data.frame")
	##df$range <- rep(i, nrow(df))##mm.df$query
	##dfList[[i]] <- df
	##df$range <- range.index
	df$range <- paste("[", range.index, "] chr ", chromosome(range), ", ID: ", sampleNames(range), sep="")
	return(df)
}
