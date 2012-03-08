.xyplot2 <- function(x, data, range, frame=50e3L, ...){
	df <- foreach(i=seq_len(nrow(range)), .combine="rbind") %do% {
		dataFrameFromRange(range=range[i, ], object=data,
				   frame=frame, range.index=i)
	}
	df$range <- factor(df$range, ordered=TRUE, levels=unique(df$range))
	df$id <- factor(df$id, ordered=TRUE, levels=unique(df$id))
	if("return.data.frame" %in% names(list(...))){
		return.df <- list(...)[["return.data.frame"]]
		if(return.df) return(df)
	}
	list.x <- as.character(x)
	i <- grep("|", list.x, fixed=TRUE)
	if(length(i) > 0){
		zvar <- list.x[[i]]
		zvar <- strsplit(zvar, " | ", fixed=T)[[1]][[2]]
		if(zvar == "range"){
			tmp <- tryCatch(df$range <- mm.df$query, error=function(e) NULL)
		}
	}
	if("gt" %in% colnames(df)){
		xyplot(x, df,
		       range=range,
		       id=df$id,
		       gt=df$gt,
		       is.snp=df$is.snp,
		       ...)
	} else {
		xyplot(x, df,
		       id=df$id,
		       range=range,
		       is.snp=df$is.snp,
		       ...)
	}
}

setMethod("xyplot2", signature(x="formula",
			       data="gSet",
			       range="RangedDataCNV"),
	  function(x, data, range, frame=50e3L, ...){
		  .xyplot2(x=x, data=data, range=range, frame=frame, ...)
	  })

setMethod("xyplot2", signature(x="formula",
			       data="SnpSet",
			       range="RangedDataCNV"),
	  function(x, data, range, frame=50e3L, ...){
		  .xyplot2(x=x, data=data, range=range, frame=frame, ...)
	  })

setMethod("xyplot", signature(x="formula", data="BeadStudioSet"),
	  function(x, data, ...){
		  if("range" %in% names(list(...))){
			  xyplot2(x, data, ...)
		  } else {
			  callNextMethod()
		  }
})

setMethod("xyplot2", signature(x="formula", data="CNSet", range="RangedDataCNV"),
	  function(x, data, range, frame=50e3L, ...){
		  z <- findOverlaps(range, data, maxgap=frame)
		  mm <- as.matrix(z)
		  mm.df <- data.frame(mm)
		  mm.df$featureNames <- featureNames(data)[mm.df$subject]
		  marker.index <- unique(mm.df$subject)
		  ##marker.index <- featuresInRange(data, rd, FRAME=frame)
		  sample.index <- match(unique(sampleNames(range)), sampleNames(data))
		  data <- data[marker.index, sample.index]
		  ## now we need to know the indices of
		  ## each range after subsetting
		  ## mm.df$subject <- match(mm.df$featureNames, featureNames(data))
		  ## we assume that each range is from a different sample
		  oligoset <- as(data, "oligoSnpSet")
		  df <- as(oligoset, "data.frame")
		  df$range.index <- mm.df$query
		  xyplot(x, df,
			 range=range,
			 gt=df$gt,
			 is.snp=df$is.snp,
			 ...)
	  })

setMethod("xyplot", signature(x="formula", data="SnpSet"),
	  function(x, data, ...){
		  if("range" %in% names(list(...))){
			  xyplot2(x, data, ...)
		  } else {
			  callNextMethod()
		  }
})

xyplotLrrBaf <- function(rd, object, frame, ...){
	index <- seq_len(nrow(rd))
	df <- foreach(i=index, .combine="rbind") %do% dataFrameFromRange(range=rd[i, ],
			       object=object, frame=frame, range.index=i)
	df$range <- factor(df$range, ordered=TRUE, levels=unique(df$range))
	if("cn" %in% colnames(df)){
		xyplot(cn~x|range, data=df,
		       baf=df$baf,
		       is.snp=df$is.snp, range=rd, ...)
	} else {
		if(!"lrr" %in% colnames(df)) stop("coercion to data.frame must have a 'lrr' column or a 'cn' column, but neither are present")
		xyplot(lrr~x|range, data=df,
		       baf=df$baf,
		       is.snp=df$is.snp, range=rd, ...)
	}
}
