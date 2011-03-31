setMethod(".plotChromosome", "eSet",
	  function(object, op){
		  if(length(unique(chromosome(object))) > 1) stop(".plotChromosome should only receive one chromosome")
		  .getCytoband <- function(object, op){
			  if(op$add.cytoband){
				  ##data(cytoband)
				  pathto <- system.file("hg18", package="SNPchip")
				  cytoband <- read.table(file.path(pathto, "cytoBand.txt"), as.is=TRUE)
				  colnames(cytoband) <- c("chrom", "start", "end", "name", "gieStain")
				  cytoband <- cytoband[cytoband[, "chrom"] == paste("chr", unique(chromosome(object)), sep=""), ]
			  }  else NULL
		  }
		  cytoband <- .getCytoband(object, op)
		  .drawCytobandWrapper <- function(S, cytoband, op, j, chromosomeName){
			  if(!op$add.cytoband) return()
			  if(nrow(cytoband) > 0)
				  plotCytoband(cytoband=cytoband,
					       new=FALSE,
					       cytoband.ycoords=op$cytoband.ycoords,
					       xlim=op$xlim[chromosomeName, ],
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
		  if(op$cytoband.side == 3) cytobandOnTop <- TRUE else cytobandOnTop <- FALSE
		  if(cytobandOnTop){
			  .drawCytobandWrapper(S=ncol(object),
					       cytoband=cytoband,
					       op=op,
					       j=j,
					       chromosomeName=unique(chromosome(object)))
		  }
		  for(j in 1:ncol(object)){
			  op$main <- op$main[j]
			  ##some parameters in op may get updated in the call to .plot
			  ##e.g., the xlim
			  op <- .plot(object[, j], op=op)
			  .drawYaxis(object=object, op=op)
			  if(!is.null(op@hmmPredict)){
				  hmmPredict <- op@hmmPredict
				  op@hmmPredict <- NULL
				  x <- split(hmmPredict[, j], chromosome(hmmPredict))
				  lapply(x, plotPredictions, op=op)
				  op@hmmPredict <- hmmPredict
			  }
			  if(op$add.centromere)
				  .drawCentromere(object[, j], op)
			  if(cytobandOnTop)
				  .drawXaxis(object=object, op=op, j=j)
		  }
		  if(!cytobandOnTop){
			  .drawCytobandWrapper(S=ncol(object),
					       cytoband=cytoband,
					       op=op,
					       j=j,
					       chromosomeName=unique(chromosome(object)))
			  .drawXaxis(object=object, op=op, j=j)
		  }
		  return(op)
          })


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


.recycle <- function(x, object, missing){
	if(length(x) == nrow(object)) return(x)
	##assume using 2 colors for homozygous and hets
	if(length(x) >= 1){
		if(inherits(object, "SnpSet")){
			gt <- as.vector(calls(object))
			gt[is.na(gt)] <- 4

			if(sum(!(gt %in% 1:4)) > 0){
				warning("Changing all genotypes that are not 1, 2, 3 to the value 4")
				gt[!(gt %in% 1:4)] <- 4
			}
			##assume that colors are to be recycled
			if(length(x) > length(unique(gt))){
				##this is a bad idea
				x <- x[1:length(unique(gt))]
			}

			##number of unique calls is greater than the
			##length of the supplied colors
			if(length(x) < length(unique(gt))){
				l <- length(unique(gt)) - length(x)
				x <- c(x, rep(missing, l))
			}

			##After above steps, the number of unique
			##calls and number of supplied colors in x
			##should be the same
			if(length(unique(gt)) == length(x)){
				x <- x[sort(unique(gt))]
			} else{
				stop("length of supplied col, bg, or plotting symbols is not equal to number of unique genotype calls.")
			}
			x <- x[gt]
		} else{
			##Calls not in assay data.  Just use the first color
			x <- rep(x[1], nrow(object))
		}
	}
	x
}

setMethod("plotSnp", "eSet",
	  function(object, hmmPredict, ...){
		  if(!missing(hmmPredict)){
			  require(VanillaICE) || stop("VanillaICE package not available")
			  i <- match(featureNames(object), featureNames(hmmPredict))
			  j <- match(sampleNames(object), sampleNames(hmmPredict))
			  hmmPredict <- hmmPredict[i, j]
		  }
		  ## create an appropriate class according to the class of
		  ## eSet
		  gp <- switch(class(object),
			       oligoSnpSet=new("ParSnpSet", snpset=object, ...),
			       SnpCallSet=new("ParSnpCallSet", snpset=object, ...),
			       SnpCopyNumberSet=new("ParSnpCopyNumberSet", snpset=object, ...),
			       RatioSnpSet=new("ParSnpSet", snpset=object, ...),
			       stop("Object is not one of the available classes"))
		  if(!missing(hmmPredict)) gp@hmmPredict <- hmmPredict
		  gp <- getPar(gp)
		  return(gp)
	  })

setMethod("plot", "eSet",
	  function(x, y, ...){
		  if(!missing(y)){
			  require(VanillaICE) || stop("VanillaICE package not available")
			  i <- match(featureNames(x), featureNames(y))
			  j <- match(sampleNames(x), sampleNames(y))
			  y <- y[i, j]
		  }
		  ## create an appropriate class according to the class of
		  ## eSet
		  gp <- switch(class(x),
			       oligoSnpSet=new("ParSnpSet", snpset=x, ...),
			       SnpCallSet=new("ParSnpCallSet", snpset=x, ...),
			       SnpCopyNumberSet=new("ParSnpCopyNumberSet", snpset=x, ...),
			       RatioSnpSet=new("ParSnpSet", snpset=x, ...),
			       stop("Object is not one of the available classes"))
		  if(!missing(y)) gp@hmmPredict <- y
		  gp <- getPar(gp)
		  return(gp)
	  })



##could we extend the trellis class? (probably too complicated)
setMethod("show", "ParESet",
	  function(object){
		  snpset <- object@snpset
		  snpset <- snpset[!is.na(chromosome(snpset)), ]
		  if(object$useLayout){
			  layout(mat=object$mat,
				 widths=object$widths,
				 heights=object$heights,
				 respect=object$respect)
		  }
		  snpList <- split(snpset, chromosome(snpset))
		  names(snpList)[names(snpList) == "X"] <- "23"
		  names(snpList)[names(snpList) == "XY"] <- "24"
		  names(snpList)[names(snpList) == "Y"] <- "25"
		  names(snpList)[names(snpList) == "M"] <- "26"
		  snpList <- snpList[order(as.numeric(names(snpList)))]
		  names(snpList)[names(snpList) == "23"] <- "X"
		  names(snpList)[names(snpList) == "24"] <- "XY"
		  names(snpList)[names(snpList) == "25"] <- "Y"
		  names(snpList)[names(snpList) == "26"] <- "M"
		  if(length(snpList) > 10) object$abbreviateChromosomeNames <- TRUE else object$abbreviateChromosomeNames <- FALSE
		  if(object$ylab == "copy number"){
			  if(any(apply(copyNumber(snpset), 2, "median", na.rm=TRUE) > 3) | any(apply(copyNumber(snpset), 2, "median", na.rm=TRUE) < 0)){
				  warning("The default ylabel 'copy number' may not be consistent with the quantity plotted on the vertical axes.  Typically, the median copy number is approximately 2 for autosomes or 1 for the male chromosome X")
			  }
		  }
		  par(allPlots(object))
		  for(i in 1:length(snpList)){
			  if(i == 1) par(yaxt="s") else par(yaxt="n")
			  object <- .plotChromosome(snpList[[i]], op=object)
		  }
		  if(object$outer.ylab) mtext(object$ylab, side=object$side.ylab, outer=TRUE, las=3, cex=object$cex.ylab, line=object$line.ylab)
		  mtext(object$xlab, side=object$side.xlab, outer=object$outer.xlab, cex=object$cex.xlab, line=object$line.xlab, adj=0)
		  if(length(sampleNames(snpset) == 1)){
			  mtext(object$main, side=3, outer=TRUE, cex=1.4)
		  }
		  object@snpset <- snpset
		  return(object)
          })


setMethod("show", "ParSnpCopyNumberSet",
	  function(object){
		  old.par <- par(no.readonly=TRUE)
		  on.exit(old.par)
		  callNextMethod()
	  })

setMethod("show", "ParSnpCallSet",
          function(object){
		  old.par <- par(no.readonly=TRUE)
		  on.exit(old.par)
		  callNextMethod()
          })

setMethod("show", "ParSnpSet",
          function(object){
		  old.par <- par(no.readonly=TRUE)
		  on.exit(old.par)
		  callNextMethod()
          })

.calculateYlim <- function(object, op){
	if("copyNumber" %in% ls(assayData(object)) | "ratio" %in% ls(assayData(object))){
		##only print this if there is more than one sample or more than 1 chromosome to plot
		if(!op$one.ylim){
			if(length(unique(chromosome(object))) > 1 || ncol(object) > 1){
				print("one.ylim is FALSE. Each panel has different ylim")
			}
		} else{
			ylim <- range(copyNumber(object), na.rm=TRUE)
		}
	} else {
		##ylimits for genotypes??
		y <- .getY(object)
		##jitter the genotype calls
		y <- jitter(y, amount=0.05)
		ylim <- range(y)
	}
	ylim
}

.plot <- function(object, op, ...){
	##could use a switch in the following statement to generate a
	##class-specific par object
	if(missing(op)) op <- new("ParESet")
	chromosomeName <- as.character(unique(chromosome(object)))

	## this changes the arguments in graphical parameters.  Should move
	## this to plotSnp().  op$ylim should be a named list
	if(!op$one.ylim){
		op$ylim <- .calculateYlim(object)
	}

	if(op$use.chromosome.size){
		op$xlim[chromosomeName, ] <- c(0, chromosomeSize(chromosomeName))
	}

	##Should do this before reordering by genotype so that the
	##returned colors are in the same order as the calls in the
	##original object
	col <- .recycle(op$col, object, missing="grey40")
	cex <- .recycle(op$cex, object, missing=op$cex[1])
	pch <- .recycle(op$pch, object, missing=op$pch[1])
	bg <- .recycle(op$bg, object, missing="grey40")
	.orderByGenotype <- function(object){
		if(!("calls" %in% ls(assayData(object)))){
			return(1:nrow(object))
		}
		gt <- as.vector(calls(object))
		gt[gt == 3] <- 1
		order(gt, decreasing=FALSE)
	}

	##ensures homozygous genotypes are plotted first
	if("calls" %in% ls(assayData(object))){
		ix <- .orderByGenotype(object)
		object <- object[ix, ]
	} else { ix <- 1:nrow(object)}
	##x <- .getX(object)##position
	x <- position(object)
	##y <- .getY(object, op)##calls or copy number
	y <- copyNumber(object)
	if(!op$outer.ylab) ylab <- op$ylab else ylab <- ""
	## this changes the arguments in graphical parameters.  Should move
	## this to plotSnp().  Might want to leave this here in case cex is
	## updated to cex=2, for instance
	##Option to recycle graphical parameters by genotype call (when available)
	plot(x=x, y=y,
	     xlim=op$xlim[chromosomeName, ],
	     ylim=op$ylim,
	     col=col[ix],
	     cex=cex[ix],
	     pch=pch[ix],
	     log=op$log,
	     bg=bg[ix],
	     xaxt="n",
	     xaxs=op$xaxs,
	     main=op$main,
	     ylab=ylab,
	     yaxt="n",
	     ...)
	if(op$abline){
		abline(h=op$abline.h, col=op$abline.col, lty=op$abline.lty, lwd=op$abline.lwd)
	}
	if(op$abline.v){
		##	  if(!is.null(op$abline.v)){
		abline(v=op$abline.v.pos, col=op$abline.v.col, lty=op$abline.v.lty, lwd=op$abline.v.lwd)
	}
	if(!is.null(op$legend)){
		legend(op$legend.location,
		       col=op$legend.col,
		       pt.bg=op$legend.bg,
		       pch=op$legend.pch,
		       legend=op$legend,
		       bty=op$legend.bty)
	}
	return(op)
}

##plotCytoband <- function(chromosome,
##                         cytoband,
##			 cytoband.ycoords,
##                         xlim,
##			 ylim=c(0, 2),
##                         xaxs="r",
##                         new=TRUE,
##                         label.cytoband=TRUE,  ##whether to label cytobands
##			 label.y=NULL,         ##if specified, use text() rather than axis()
##			 srt,
##                         cex.axis=1,
##                         outer=FALSE,
##			 taper=0.15,
##                         ...){
##	def.par <- par(no.readonly=TRUE)
##	on.exit(def.par)
##	if(missing(cytoband)) data(cytoband, package="SNPchip", envir=environment())
##	if(missing(chromosome)){
##		if(length(unique(cytoband[, "chrom"])) > 1) stop("Must specify chromosome")
##	}
##	if(length(unique(cytoband$chrom)) > 1){
##		cytoband <- cytoband[cytoband[, "chrom"] == chromosome, ]
##	}
##	if(missing(cytoband.ycoords)){
##		cytoband.ycoords <- ylim
##	}
##	rownames(cytoband) <- as.character(cytoband[, "name"])
##	if(missing(xlim)) xlim <- c(0, chromosomeSize(unique(cytoband$chrom)))
##	cytoband_p <- cytoband[grep("^p", rownames(cytoband), value=TRUE), ]
##	cytoband_q <- cytoband[grep("^q", rownames(cytoband), value=TRUE), ]
##
##	p.bands <- nrow(cytoband_p)
##	cut.left  <- c()
##	cut.right <- c()
##	##  1st  band of arm or 1st  band after  "stalk"
##	##  last band of arm or last band before "stalk"
##	for (i in 1:nrow(cytoband)) {
##		if (i == 1)                             { cut.left[i] <- TRUE; cut.right[i] <- FALSE} else
##		if (i == p.bands)                       { cut.left[i] <- FALSE; cut.right[i] <- TRUE} else
##		if (i == (p.bands+1))                   { cut.left[i] <- TRUE; cut.right[i] <- FALSE} else
##		if (i == nrow(cytoband)) { cut.left[i] <- FALSE; cut.right[i] <- TRUE} else{
##			cut.left[i] <- FALSE; cut.right[i] <- FALSE
##		}
##	}
##	for (i in 1:nrow(cytoband)) {
##		if (as.character(cytoband[i, "gieStain"]) == "stalk") {
##			cut.right[i-1] <- TRUE
##			cut.left[i] <- NA
##			cut.right[i] <- NA
##			cut.left[i+1] <- TRUE
##		}
##	}
##	##When plotting subregions of a chromosome, this prevents the
##	##cytobands from extending beyond the subsetted object
##
##	## exclude cytobands that end before the minimum plotting
##	##limits
##	include <- cytoband[, "chromEnd"] > xlim[1] & cytoband[, "chromStart"] < xlim[2]
##
##	cytoband <- cytoband[include, ]
##	N <- nrow(cytoband)
##	cytoband[N, "chromEnd"] <- min(xlim[2], cytoband[N, "chromEnd"])
##	cytoband[1, "chromStart"] <- max(xlim[1], cytoband[1, "chromStart"])
##	cut.left <- cut.left[include]
##	cut.right <- cut.right[include]
##	if(new){
##		xx <- c(0, cytoband[nrow(cytoband), "chromEnd"])
##		yy <- c(0, 2)
####		yy <- ylim
##		plot(xx,
##		     yy,
##		     xlim=xlim,
##		     type="n",
##		     xlab="",
##		     ylab="",
##		     axes=FALSE,
##		     xaxs=xaxs,
##		     ...)
##	}
##	top <- cytoband.ycoords[2]
##	bot <- cytoband.ycoords[1]
##	h <- top-bot
##	p <- taper
##	for (i in 1:nrow(cytoband)) {
##		start <- cytoband[i, "chromStart"]
##		last   <- cytoband[i, "chromEnd"]
##		delta = (last-start)/4
##		getStain <- function(stain){
##			switch(stain,
##			       gneg="grey100",
##			       gpos25="grey90",
##			       gpos50="grey70",
##			       gpos75="grey40",
##			       gpos100="grey0",
##			       gvar="grey100",
##			       stalk="brown3",
##			       acen="brown4",
##			       "white")
##		}
##		color <- getStain(as.character(cytoband[i, "gieStain"]))
##		if (is.na(cut.left[i]) & is.na(cut.right[i])) {
##			## this is a "stalk", do not draw box. Draw two vertical lines instead
##			delta <- (last-start)/3
##			lines(c(start+delta, start+delta), ylim, col=color)
##			lines(c(last-delta, last-delta), ylim, col=color)
##		} else if (cut.left[i] & cut.right[i]) {      # cut both lasts
####			polygon(c(start, start+delta, last-delta, last, last, last-delta, start+delta, start),
####				c(0.3, 0, 0, 0.3, 1.7, 2, 2, 1.7), col=color)
##			##Taper both ends
##			yy <- c(bot + p*h, bot, bot, bot + p*h, top - p*h, top, top, top - p*h)
##			polygon(c(start, start+delta, last-delta, last, last, last-delta, start+delta, start),
##				yy, col=color)
##		} else if (cut.left[i]) {              # cut left last only
##			##Taper left end only
##			yy <- c(bot + p*h, bot, bot, top, top, top - p*h)
##			polygon(c(start, start+delta, last, last, start+delta, start),
##				yy, col=color)
##		} else if (cut.right[i]) {             # cut right last only
##			##Taper right end only
##			yy <- c(bot, bot, bot + p*h, top - p*h, top, top)
##			polygon(c(start, last-delta, last, last, last-delta, start),
##				yy,col=color)
##		} else {
##			##Rectangle
##			polygon(c(start, last, last, start),
##				c(bot, bot, top, top), col=color)
##		}
##	}
##	my.x <- (cytoband$chromStart+cytoband$chromEnd)/2
##	if(label.cytoband){
##		if(is.null(label.y)){
##			##if plotting on a new device
##			axis(1,
##			     at=my.x,
##			     labels=rownames(cytoband),
##			     outer=outer,
##			     cex.axis=cex.axis,
##			     line=1,
##			     las=3)
##		} else{
##			##put cytoband labels at height label.y
##			if(!is.numeric(label.y)){
##				warning("label.y must be numeric -- using default y coordinates for cytoband labels")
##				label.y <- bot - p*h
##			}
##			if(missing(srt)) srt <- 90
##			text(x=my.x,
##			     y=rep(label.y, length(my.x)),
##			     labels=rownames(cytoband),
##			     srt=srt)
##		}
##	}
##	return()
##}

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


plotPredictions <- function(object, op){
	breakpoints <- breakpoints(object)
	position <- position(object)
	chromosome <- chromosome(object)
	states <- states(object)
	predictions <- predictions(object)

	chr <- unique(chromosome)
	for(i in chr){
	  tmp <- breakpoints[breakpoints[, "chr"] == i, , drop=FALSE]
	  .drawRect <- function(x, position, op){
		  col <- op$col.predict
		  ##if(length(x) < 1) return()
		  start <- max(as.numeric(x["start"]), op$xlim[1])
		  last <- min(as.numeric(x["end"]), op$xlim[2])
		  predict <- predictions[position >= start & position <= last]
		  predict <- predict[!is.na(predict)]
		  if(length(unique(predict)) > 1) {
			  stop("predictions not unique")
		  }
		  col <- col[unique(predict)]
		  rect(xleft=start,
		       ybottom=op$hmm.ycoords[1],
		       xright=last,
		       ytop=op$hmm.ycoords[2],
		       col=col,
		       border=col)
	  }
	  if(!is.null(tmp)){
		  ##best to draw the biggest regions first and the smallest regions last.
		  tmp <- tmp[order(tmp[, "end"]-tmp[, "start"], decreasing=TRUE), ]
		  apply(tmp, 1, .drawRect, position=position, op=op)
	  }
	  if(!is.null(op$legend.predict)){
		  legend(op$legend.location.predict,
			 fill=op$col.predict,
			 legend=op$legend.predict,
			 ncol=length(op$col.predict),
			 bty=op$legend.bty,
			 cex=op$cex.legend)
	  }
  }
	return(op)
}
