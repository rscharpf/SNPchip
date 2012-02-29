centromere <- function(chromosome, build="hg18", verbose=FALSE){
	if(verbose)  message(paste("centromere coordinates based on build", build))
	if(missing(chromosome) | !(chromosome %in% c(1:22, "X"))) stop("must specify chromosome 1-22, or X as character string")
	pathto <- system.file(build, package="SNPchip")
	tmp <- read.table(file.path(pathto, "centromeres.txt"), as.is=TRUE)
	as.integer(tmp[paste("chr", chromosome, sep=""), ])
}


chromosomeSize <- function(chromosome, build="hg18", verbose=FALSE){
	if(verbose) message(paste("chromosome size using build", build))
	if(!is.character(chromosome)) stop("argument to chromosomeSize must be one of the following character strings: 1, ..., 22, X, or Y")
	if(length(grep("chr", chromosome)) == 0) chromosome <- paste("chr", chromosome, sep="")
	if(any(!(chromosome %in% paste("chr", c(1:22, "X", "Y", "XY", "M"), sep="")))) stop("chromosome must be chr1-22, chrX, chrY, or chrM")
	pathto <- system.file(build, package="SNPchip")
	tmp <- read.table(file.path(pathto, "chromInfo.txt"), as.is=TRUE, row.names=1)
	##data(chromosomeAnnotation, package="SNPchip", envir=environment())
	tmp[chromosome, 1]
}

.getCytoband <- function(object, op){
	##browser()
	if(op$add.cytoband){
		##data(cytoband)
		pathto <- system.file("hg18", package="SNPchip")
		cytoband <- read.table(file.path(pathto, "cytoBand.txt"), as.is=TRUE)
		colnames(cytoband) <- c("chrom", "start", "end", "name", "gieStain")
		cytoband <- cytoband[cytoband[, "chrom"] == paste("chr", unique(chromosome(object)), sep=""), ]
	}  else NULL
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

plotCytoband <- function(chromosome,
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
			 build="hg18",
                         ...){
	##def.par <- par(no.readonly=TRUE)
	##on.exit(def.par)
	if("use.lattice" %in% names(list(...))) {
		segments <- lsegments
		polygon <- lpolygon
	}
	if(missing(cytoband)){
		if(verbose) message(paste("Cytoband annotation obtained from build", build))
		pathto <- system.file("hg18", package="SNPchip")
		cytoband <- read.table(file.path(pathto, "cytoBand.txt"), as.is=TRUE)
		colnames(cytoband) <- c("chrom", "start", "end", "name", "gieStain")
		##data(cytoband, package="SNPchip", envir=environment())
	}
	if(missing(chromosome)){
		if(length(unique(cytoband[, "chrom"])) > 1) stop("Must specify chromosome")
	}
	if(length(unique(cytoband$chrom)) > 1){
		cytoband <- cytoband[cytoband[, "chrom"] == paste("chr", chromosome, sep=""), ]
	}
	if(missing(cytoband.ycoords)){
		cytoband.ycoords <- ylim
	}
	rownames(cytoband) <- as.character(cytoband[, "name"])
	if(missing(xlim)) xlim <- c(0, chromosomeSize(unique(cytoband$chrom)))
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
	if(label.cytoband){
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

