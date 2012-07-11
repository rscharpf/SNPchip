centromere <- function(chromosome, build, verbose=FALSE){
	chr <- cleanname(chromosome)
	if(missing(build)) stop("UCSC genome build must be specified (hg18 or hg19 supported)")
	if(!build %in% c("hg18", "hg19")) stop("Only hg18 and hg19 UCSC builds supported")
	if(verbose)  message(paste("centromere coordinates based on build", build))
	pathto <- system.file("extdata", package="SNPchip")
	gaps <- readRDS(file.path(pathto, paste("gap_", build, ".rda", sep="")))
	if(!chr %in% chromosome(gaps)) stop(paste("arg chromosome must be one of ", chromosome(gaps)), sep="")
	gap <- gaps[match(chr, chromosome(gaps)), ]
	c(start(gap), end(gap))
}


chromosomeSize <- function(chromosome, build, verbose=FALSE){
	if(missing(build)) stop("UCSC genome build must be specified (hg18 or hg19 supported)")
	chromosome <- cleanname(chromosome)
	getSequenceLengths(build)[chromosome]
}



.drawCytobandWrapper <- function(S, cytoband, op, j, chromosomeName){
	if(!op$add.cytoband) return()
	if(nrow(cytoband) > 0)
		plotCytoband(cytoband=cytoband,
			     new=FALSE,
			     cytoband.ycoords=op$cytoband.ycoords,
			     xlim=op$xlim[as.character(chromosomeName), ],
			     xaxs=op$xaxs,
			     label.cytoband=op$label.cytoband,
			     srt=op$cytoband.srt,
			     label.y=op$cytoband.label.y,
			     cex.axis=op$cex.axis,
			     outer=op$outer.cytoband.axis,
			     taper=op$cytoband.taper)
}

.drawYaxis <- function(object, op, j){
	if(unique(chromosome(object)) != op$firstChromosome) return()
	if(op$yaxt == "n") return()
	if("copyNumber" %in% ls(assayData(object)) | "ratio" %in% ls(assayData(object))){
		at <- pretty(op$ylim)
		at <- c(op$at, at)
		labels <- at
	}  else {
		at <- c(0, 1)
		labels <- c("AA/BB", "AB")
	}
	axis(side=2, at=at, labels=labels, las=1, cex.axis=op$cex.axis)
}

.drawCentromere <- function(object, op){
	centromere.coords <- centromere(unique(chromosome(object)))
	##data(chromosomeAnnotation, package="SNPchip", envir=environment())
	##centromere <- chromosomeAnnotation[unique(chromosome(object)), ]
	xleft <- centromere.coords[1]
	xright <- centromere.coords[2]
	rect(xleft=xleft, ybottom=op$ylim[1],
	     xright=xright, ytop=op$ylim[2],
	     col=op$col.centromere,
	     border=op$border.centromere)
}



.drawXaxis <- function(object, op, j){
	chromosomeName <- as.character(unique(chromosome(object)))
	if(op$xaxt == "n") return()
	if(op$alternate.xaxis.side){
		side <- op$xaxis.side[[unique(chromosome(object))]]
	} else side <- op$xaxis.side
	if(side == 1 & j == ncol(object) | side == 3 & j == 1){
		##if labeling the cytoband, force the x-axis to be drawn on top
		if(op$cytoband.side==1 & op$label.cytoband){
##			op$outer.axis <- TRUE
##			op$line.axis <- 2
			side <- 3
		}
		axis(side,
		     at=pretty(op$xlim[chromosomeName, ], op$lab[2]),
		     outer=op$outer.axis,
		     labels=pretty(op$xlim[chromosomeName, ], op$lab[2])/1e6,
		     cex.axis=op$cex.axis,
		     col=op$col.axis,
		     col.axis=op$col.axis,
		     las=1,
		     line=op$line.axis,
		     lwd=1,
		     mgp=c(2, 0.5, 0))
		if(op$label.chromosome){
			if(!op$abbreviateChromosomeNames){
				chr <- paste("Chr", unique(chromosome(object)))
			} else{
				chr <- unique(chromosome(object))
			}
			mtext(chr,
			      side=side,
			      outer=FALSE,
			      line=op$line.label.chromosome,
			      cex=op$cex.lab)
		}
	}
}

cleanname <- function(chromosome){
	if(is.numeric(chromosome)) {
		chromosome <- paste("chr",integer2chromosome(chromosome), sep="")
	} else {
		x <- strsplit(chromosome, "chr")[[1]]
		if(length(x)==1) chromosome <- paste("chr", x, sep="")
	}
	return(chromosome)
}

plotIdiogram <- function(chromosome,
			 build,
                         cytoband,
			 cytoband.ycoords,
                         xlim,
			 ylim=c(0, 2),
                         new=TRUE,
                         label.cytoband=TRUE,  ##whether to label cytobands
			 label.y=NULL,         ##if specified, use text() rather than axis()
			 srt,
                         cex.axis=1,
                         outer=FALSE,
			 taper=0.15,
			 verbose=FALSE,
			 unit=c("bp", "Mb"),
			 is.lattice=FALSE,
                         ...){
	##def.par <- par(no.readonly=TRUE)
	##on.exit(def.par)
	if(missing(build)) stop("must specify genome build")
	if(is.lattice){
		segments <- lsegments
		polygon <- lpolygon
	}
	if(missing(cytoband)){
		pathto <- system.file("extdata", package="SNPchip")
		if(verbose) message("Reading cytoband annotation for UCSC genome build ", build)
		cytoband <- read.table(file.path(pathto, paste("cytoBand_", build, ".txt", sep="")), as.is=TRUE)
		colnames(cytoband) <- c("chrom", "start", "end", "name", "gieStain")
	}
	if(!missing(chromosome)){
		chromosome <- cleanname(chromosome)
	} else {
		if(length(unique(cytoband[, "chrom"])) > 1) stop("Must specify chromosome")
	}
	##if(length(unique(cytoband$chrom)) > 1){
	cytoband <- cytoband[cytoband[, "chrom"] == chromosome, ]
	unit <- match.arg(unit)
	if(unit=="Mb"){
		cytoband$start <- cytoband$start/1e6
		cytoband$end <- cytoband$end/1e6
	}
	if(missing(cytoband.ycoords)){
		cytoband.ycoords <- ylim
	}
	rownames(cytoband) <- as.character(cytoband[, "name"])
	sl <- getSequenceLengths(build)[chromosome]
	if(missing(xlim)) xlim <- c(0, sl)
	if(unit=="Mb") xlim <- xlim/1e6
	cytoband_p <- cytoband[grep("^p", rownames(cytoband), value=TRUE), ]
	cytoband_q <- cytoband[grep("^q", rownames(cytoband), value=TRUE), ]

	p.bands <- nrow(cytoband_p)
	cut.left  <- c()
	cut.right <- c()
	##  1st  band of arm or 1st  band after  "stalk"
	##  last band of arm or last band before "stalk"
	for (i in 1:nrow(cytoband)) {
		if (i == 1)                             { cut.left[i] <- TRUE; cut.right[i] <- FALSE} else
		if (i == p.bands)                       { cut.left[i] <- FALSE; cut.right[i] <- TRUE} else
		if (i == (p.bands+1))                   { cut.left[i] <- TRUE; cut.right[i] <- FALSE} else
		if (i == nrow(cytoband)) { cut.left[i] <- FALSE; cut.right[i] <- TRUE} else{
			cut.left[i] <- FALSE; cut.right[i] <- FALSE
		}
	}
	for (i in 1:nrow(cytoband)) {
		if (as.character(cytoband[i, "gieStain"]) == "stalk") {
			cut.right[i-1] <- TRUE
			cut.left[i] <- NA
			cut.right[i] <- NA
			cut.left[i+1] <- TRUE
		}
	}
	##When plotting subregions of a chromosome, this prevents the
	##cytobands from extending beyond the subsetted object
	##
	## exclude cytobands that end before the minimum plotting
	##limits
	include <- cytoband[, "end"] > xlim[1] & cytoband[, "start"] < xlim[2]
	cytoband <- cytoband[include, ]
	N <- nrow(cytoband)
	cytoband[N, "end"] <- min(xlim[2], cytoband[N, "end"])
	cytoband[1, "start"] <- max(xlim[1], cytoband[1, "start"])
	cut.left <- cut.left[include]
	cut.right <- cut.right[include]
	if(new){
		xx <- c(0, cytoband[nrow(cytoband), "end"])
		yy <- cytoband.ycoords
##		yy <- ylim
		plot(xx,
		     yy,
		     xlim=xlim,
		     type="n",
		     xlab="",
		     ylab="",
		     axes=FALSE,
		     yaxs="i",
		     ylim=ylim,
		     ...)
	}
	top <- cytoband.ycoords[2]
	bot <- cytoband.ycoords[1]
	h <- top-bot
	p <- taper
	for(i in seq_len(nrow(cytoband))) {
		start <- cytoband[i, "start"]
		last   <- cytoband[i, "end"]
		delta = (last-start)/4
		getStain <- function(stain){
			switch(stain,
			       gneg="grey100",
			       gpos25="grey90",
			       gpos50="grey70",
			       gpos75="grey40",
			       gpos100="grey0",
			       gvar="grey100",
			       stalk="brown3",
			       acen="brown4",
			       "white")
		}
		color <- getStain(as.character(cytoband[i, "gieStain"]))
		if (is.na(cut.left[i]) & is.na(cut.right[i])) {
			## this is a "stalk", do not draw box. Draw two vertical lines instead
			delta <- (last-start)/3
			segments(start+delta, cytoband.ycoords[1], start+delta, cytoband.ycoords[2])
			segments(last-delta, cytoband.ycoords[1], last-delta, cytoband.ycoords[2])
##			lines(c(start+delta, start+delta), ylim, col=color)
##			lines(c(last-delta, last-delta), ylim, col=color)
		} else if (cut.left[i] & cut.right[i]) {      # cut both lasts
##			polygon(c(start, start+delta, last-delta, last, last, last-delta, start+delta, start),
##				c(0.3, 0, 0, 0.3, 1.7, 2, 2, 1.7), col=color)
			##Taper both ends
			yy <- c(bot + p*h, bot, bot, bot + p*h, top - p*h, top, top, top - p*h)
			polygon(c(start, start+delta, last-delta, last, last, last-delta, start+delta, start),
				yy, col=color)
		} else if (cut.left[i]) {              # cut left last only
			##Taper left end only
			yy <- c(bot + p*h, bot, bot, top, top, top - p*h)
			polygon(c(start, start+delta, last, last, start+delta, start),
				yy, col=color)
		} else if (cut.right[i]) {             # cut right last only
			##Taper right end only
			yy <- c(bot, bot, bot + p*h, top - p*h, top, top)
			polygon(c(start, last-delta, last, last, last-delta, start),
				yy,col=color)
		} else {
			##Rectangle
			polygon(c(start, last, last, start),
				c(bot, bot, top, top), col=color)
		}
	}
	my.x <- (cytoband[, "start"] + cytoband[, "end"])/2
	if(label.cytoband & !is.lattice){
		if(is.null(label.y)){
			##if plotting on a new device
			axis(1,
			     at=my.x,
			     labels=rownames(cytoband),
			     outer=outer,
			     cex.axis=cex.axis,
			     line=1,
			     las=3)
		} else{
			##put cytoband labels at height label.y
			if(!is.numeric(label.y)){
				warning("label.y must be numeric -- using default y coordinates for cytoband labels")
				label.y <- bot - p*h
			}
			if(missing(srt)) srt <- 90
			text(x=my.x,
			     y=rep(label.y, length(my.x)),
			     labels=rownames(cytoband),
			     srt=srt)
		}
	}
	return()
}

plotCytoband <- function(...) .Deprecated("plotCytoband is deprecated. Use plotIdiogram instead.")

plotCytoband2 <- function(chromosome,
			  build,
			  cytoband,
			  xlim,
			  xaxs="r",
			  new=TRUE,
			  label.cytoband=TRUE,
			  cex.axis=1,
			  outer=FALSE,
			  verbose=TRUE,
			  ...){
	def.par <- par(no.readonly=TRUE, mar=c(4.1, 0.1, 3.1, 2.1))
	on.exit(def.par)
	if(missing(build)) stop("must specify genome build")
##	if(is.lattice){
##		segments <- lsegments
##		polygon <- lpolygon
##	}
	if(missing(cytoband)){
		pathto <- system.file("extdata", package="SNPchip")
		if(verbose) message("Reading cytoband annotation for UCSC genome build ", build)
		cytoband <- read.table(file.path(pathto, paste("cytoBand_", build, ".txt", sep="")), as.is=TRUE)
		colnames(cytoband) <- c("chrom", "start", "end", "name", "gieStain")
	}
	if(!missing(chromosome)){
		chromosome <- cleanname(chromosome)
	} else stop("must specify chromosome")
	cytoband <- cytoband[cytoband[, "chrom"] == chromosome, ]
	rownames(cytoband) <- as.character(cytoband[, "name"])
	sl <- getSequenceLengths(build)[chromosome]
	if(missing(xlim)) xlim <- c(0, sl)
	cytoband_p <- cytoband[grep("^p", rownames(cytoband), value=TRUE), ]
	cytoband_q <- cytoband[grep("^q", rownames(cytoband), value=TRUE), ]

	p.bands <- nrow(cytoband_p)
	cut.left  <- c()
	cut.right <- c()
	##  1st  band of arm or 1st  band after  "stalk"
	##  last band of arm or last band before "stalk"
	for (i in 1:nrow(cytoband)) {
		if (i == 1) { cut.left[i] <- TRUE; cut.right[i] <- FALSE} else
		if (i == p.bands) { cut.left[i] <- FALSE; cut.right[i] <- TRUE} else
		if (i == (p.bands+1)) { cut.left[i] <- TRUE; cut.right[i] <- FALSE} else
		if (i == nrow(cytoband)) { cut.left[i] <- FALSE; cut.right[i] <- TRUE} else{
			cut.left[i] <- FALSE; cut.right[i] <- FALSE
		}
	}
	for (i in 1:nrow(cytoband)) {
		if (as.character(cytoband[i, "gieStain"]) == "stalk") {
			cut.right[i-1] <- TRUE
			cut.left[i] <- NA
			cut.right[i] <- NA
			cut.left[i+1] <- TRUE
		}
	}
	##When plotting subregions of a chromosome, this prevents the cytobands from extending beyond the subsetted object
	##exclude cytobands that end before the minimum plotting limits
	include <- cytoband[, "end"] > xlim[1] & cytoband[, "start"] < xlim[2]
	cytoband <- cytoband[include, ]
	cut.left <- cut.left[include]
	cut.right <- cut.right[include]
	if(new){
		plot(c(0, cytoband[nrow(cytoband), "end"]),
		     c(0, 2),
		     xlim=xlim,
		     type="n",
		     xlab="",
		     ylab="",
		     axes=FALSE,
		     xaxs=xaxs,
		     ...)
	}
	for (i in 1:nrow(cytoband)) {
		start <- cytoband[i, "start"]
		last   <- cytoband[i, "end"]
		delta = (last-start)/4
		getStain <- function(stain){
			switch(stain,
			       gneg="grey100",
			       gpos25="grey90",
			       gpos50="grey70",
			       gpos75="grey40",
			       gpos100="grey0",
			       gvar="grey100",
			       stalk="brown3",
			       acen="brown4",
			       "white")
		}
		color <- getStain(as.character(cytoband[i, "gieStain"]))
		if (is.na(cut.left[i]) & is.na(cut.right[i])) {
			## this is a "stalk", do not draw box. Draw two vertical lines instead
			delta <- (last-start)/3
			lines(c(start+delta, start+delta), c(0,2), col=color)
			lines(c(last-delta, last-delta), c(0,2), col=color)
		} else if (cut.left[i] & cut.right[i]) {      # cut both lasts
			polygon(c(start, start+delta, last-delta, last, last, last-delta, start+delta, start),
				c(0.3, 0, 0, 0.3, 1.7, 2, 2, 1.7), col=color)
		} else if (cut.left[i]) {              # cut left last only
			polygon(c(start, start+delta, last, last, start+delta, start),
				c(0.3, 0, 0, 2, 2, 1.7), col=color)
		} else if (cut.right[i]) {             # cut right last only
			polygon(c(start, last-delta, last, last, last-delta, start),
				c(0, 0, 0.3, 1.7, 2, 2),col=color)
		} else {
			polygon(c(start, last, last, start),
				c(0, 0, 2, 2), col=color)
		}
	}
	my.x <- (cytoband$start+cytoband$end)/2
	if(label.cytoband){
		axis(1, at=my.x,
		     labels=rownames(cytoband),
		     outer=outer,
		     cex.axis=cex.axis,
		     line=1, las=3, tick=FALSE)
		axis(1, at=cytoband$start,
		     outer=outer,
		     cex.axis=cex.axis,
		     line=1, las=3, label=FALSE)
	}
	return()
}

getCytoband <- function(build){
	path <- system.file("extdata", package="SNPchip")
	if (missing(build)) build <- "hg19"
	cytoband <- read.table(file.path(path, paste("cytoBand_", build, ".txt", sep="")), as.is=TRUE, header=FALSE)
	colnames(cytoband) <- c("chrom", "start", "end", "name", "gieStain")
	return(cytoband)
}
