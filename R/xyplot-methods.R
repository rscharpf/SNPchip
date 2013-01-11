.xyplot2 <- function(x, data, range, frame=50e3L, panel, ...){
	L <- length(start(range))
	df <- foreach(i=seq_len(L), .combine="rbind") %do% {
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
##		if(zvar == "range"){
##			tmp <- tryCatch(df$range <- mm.df$query, error=function(e) NULL)
##		}
	}
	if("gt" %in% colnames(df)){
		lattice:::xyplot(x, df,
				range=range,
				id=df$id,
				gt=df$gt,
				is.snp=df$is.snp,
				panel=panel,
				...)
	} else {
		lattice::xyplot(x, df,
				id=df$id,
				range=range,
				is.snp=df$is.snp,
				panel=panel,
				...)
	}
}

setMethod("xyplot2", signature(x="formula",
			       data="gSet"),
	  function(x, data, range, frame=50e3L, ...){
		  .xyplot2(x=x, data=data, range=range, frame=frame, ...)
	  })

setMethod("xyplot2", signature(x="formula",
			       data="SnpSet"),
	  function(x, data, range, frame=50e3L, ...){
		  .xyplot2(x=x, data=data, range=range, frame=frame, ...)
	  })

setMethod("xyplot2", signature(x="formula", data="CNSet", range="RangedDataCNV"),
	  function(x, data, range, frame=50e3L, ...){
		  if(is(range, "RangedDataCNV"))
			  z <- findOverlaps(range, data, maxgap=frame)
		  if(is(range, "GRanges")){
			  frange <- oligoClasses:::makeFeatureGRanges(object)
			  z <- findOverlaps(range, frange)
		  }
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

setMethod("xyplot", signature(x="formula", data="BeadStudioSet"),
	  function(x, data, ...){
		  if("range" %in% names(list(...))){
			  xyplot2(x, data, ...)
		  } else {
			  callNextMethod()
		  }
})

setMethod("xyplot", signature(x="formula", data="gSet"),
	  function(x, data, ...){
		  if("range" %in% names(list(...))){
			  ##if(is(list(...)[["range"]], "RangedDataCNV"))
			  xyplot2(x, data, ...)
		  } else {
			  callNextMethod()
		  }
})

xyplotOligoSnpSet <- function(x, object, ...){
	df <- as(object, "data.frame")
	xyplot(x, df, ...)
}

xyplotRangeInOligoSnpSet <- function(x, object, ...){
	df <- as(object, "data.frame")
	xyplot(x, df, ...)
}

xyplotLrrBaf <- function(rd, object, frame, ...){
	if(is(rd, "GRangesList")) {
		rd <- stack(rd)
		index <- seq_len(length(rd))
	} else {
		if(is(rd, "GRanges")) index <- seq_len(length(rd))
		if(is(rd, "RangedDataCNV")) index <- seq_len(nrow(rd))
	}
	if(any(is.na(position(object)))) stop("NA values not permitted in position(object)")
	i <- NULL
	if(is(object, "gSet")){
		df <- foreach(i=index, .combine="rbind") %do% dataFrameFromRange(range=rd[i, ],
				       object=object, frame=frame, range.index=i)
	} else {
		if(!is(object, "gSetList")) stop("object must extend gSet or gSetList")
		chr <- unique(chromosome(rd))
		object <- object[paste("chr", chromosome(object), sep="") %in% chr]
		dflist <- list()
		for(i in seq_along(chr)){
			obj <- object[[i]]
			rd2 <- rd[chromosome(rd) == paste("chr", chromosome(obj)[1],sep=""),  ]
			index <- seq_len(length(rd2))
			dflist[[i]] <- foreach(i=index, .combine="rbind") %do% dataFrameFromRange(range=rd2[i, ],
							object=obj, frame=frame, range.index=i)
		}
		df <- do.call("rbind", dflist)
	}
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

##xyplotLrrBaf <- function(rd, object, frame, ...){
##	if(is(rd, "GRangesList")) {
##		rd <- stack(rd)
##		index <- seq_len(length(rd))
##	} else {
##		if(is(rd, "GRanges")) index <- seq_len(length(rd))
##		if(is(rd, "RangedDataCNV")) index <- seq_len(nrow(rd))
##	}
##	if(any(is.na(position(object)))) stop("NA values not permitted in position(object)")
##	i <- NULL
##	df <- foreach(i=index, .combine="rbind") %do% dataFrameFromRange(range=rd[i, ],
##			       object=object, frame=frame, range.index=i)
##	df$range <- factor(df$range, ordered=TRUE, levels=unique(df$range))
##	if("cn" %in% colnames(df)){
##		xyplot(cn~x|range, data=df,
##		       baf=df$baf,
##		       is.snp=df$is.snp, range=rd, ...)
##	} else {
##		if(!"lrr" %in% colnames(df)) stop("coercion to data.frame must have a 'lrr' column or a 'cn' column, but neither are present")
##		xyplot(lrr~x|range, data=df,
##		       baf=df$baf,
##		       is.snp=df$is.snp, range=rd, ...)
##	}
##}

latticeFigs <- function(gr, data, colors, ...){
	intervals <- unique(data$interval)
	lrr.fig <- xyplot(lrr~x/1e6 | interval, data=data, ...,
			  scales=list(x=list(relation="free", axs="i")),
			  ylab="log R ratios",
			  granges=gr,
			  colors=colors,
			  ylim=c(-2, 1.6),
			  id = data$id,
			  panel=function(x,y, granges, colors, id, ..., subscripts){
				  id <- id[subscripts[1]]
				  colors <- colors[sampleNames(granges) == id]
				  granges <- granges[sampleNames(granges) == id, ]
				  panel.grid(h=3,v=3)
				  panel.xyplot(x, y, ...)
				  lrect(xleft=start(granges)/1e6,
					xright=end(granges)/1e6,
					ytop=rep(1.5, length(granges)),
					ybottom=rep(1.3, length(granges)),
					col=colors,
					border=colors)
				  x1 <- current.panel.limits()[[1]]
				  ltext(x1, 1.5, "HMM track", cex=0.7, adj=0)
			  },
			  xlab="position (Mb)", layout=c(length(intervals), 1))
	baf.fig <- xyplot(baf~x/1e6 | interval, data=data, ...,
			  scales=list(x=list(relation="free", axs="i")),
			  ylab="BAFs",
			  granges=gr,
			  colors=colors,
			  ylim=c(0, 1),
			  panel=function(x,y, granges, colors, ...){
				  panel.grid(h=3, v=3)
				  panel.xyplot(x, y, ...)
			  },
			  xlab="position (Mb)", layout=c(length(intervals),1))
	figs <- list(lrr=lrr.fig, baf=baf.fig)
}

arrangeFigs <- function(lattice.figs, ...){
	## up to the panel function to plot all ranges
	##intervals <- unique(data$interval)
	lrr.fig <- lattice.figs[["lrr"]]
	baf.fig <- lattice.figs[["baf"]]
	grid.newpage()
	trellis.par.set("fontsize", list(text=10))
	lvp <- viewport(x=0,
			y=0.5,
			width=unit(1, "npc"),
			height=unit(0.5, "npc"), just=c("left", "bottom"),
			name="lvp")
	pushViewport(lvp)
	pushViewport(dataViewport(xscale=c(0,1),
				  yscale=c(0.05,1), clip="on"))
	print(lrr.fig, newpage=FALSE, prefix="plot1", more=TRUE)
	upViewport(0)
	lvp2 <- viewport(x=0, y=0,
			 width=unit(1, "npc"),
			 height=unit(0.5, "npc"),
			 just=c("left", "bottom"), name="lvp2")
	pushViewport(lvp2)
	pushViewport(dataViewport(xscale=c(0,1), yscale=c(0.05,1), clip="on"))
	print(baf.fig, newpage=FALSE, prefix="plot2", more=FALSE)
	upViewport(0)
}
