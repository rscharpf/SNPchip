##updates graphical parameters with information from the data class
setMethod("getPar", "ParESet", 
          function(object, build="hg18", ...){
		  if(!is.null(hmmPredict(object))){
			  hmmPredict <- hmmPredict(object)
			  if(is.null(object$col.predict)){
				  require(RColorBrewer, quietly=TRUE) || stop("RColorBrewer package not available")
				  if(length(states(hmmPredict)) >= 3){
					  col.predict <- brewer.pal(length(states(hmmPredict)), "BrBG")
				  } else{
					  col.predict <- rep("black", 2)
				  }
				  col.predict[states(hmmPredict) == "N"] <- "white"
				  print("col.predict not specified in list of graphical parameters. Using the following colors:")
				  print(col.predict)
				  object$col.predict <- col.predict
				  object$legend.fill.predict <- col.predict
			  }
			  if(is.null(object$height.predictions)){
				  object$height.predict <- 0.2
			  }
			  ## indicator of whether to draw vertical lines at breakpoints
			  if(object$abline.v){
				  ##if position of vertical lines are not specified
				  if(is.null(object$abline.v.pos)){
					  if(length(sampleNames(hmmPredict)) == 1){
						  v <- breakpoints(hmmPredict)[, c("state", "start", "last")]
						  v <- v[v$state != "N", ]
						  v <- c(v$start, v$last)*1e6
						  object$abline.v.pos <- v
					  }
				  }
			  }			  
		  }
		  snpset <- snpset(object)
		  snpset <- snpset[!is.na(chromosome(snpset)), ]
		  ##layout
		  chromosomeNames <- integer2chromosome(unique(chromosome(snpset)))
		  chromosomeNames <- chromosomeNames[order(chromosomeNames)]
		  N <- length(chromosomeNames)
		  S <- ncol(snpset)
		  pathto <- system.file(build, package="SNPchip")
		  chromosomeAnnotation <- read.table(file.path(pathto, "chromInfo.txt"), as.is=TRUE, row.names=1)
		  ##data(chromosomeAnnotation, package="SNPchip", envir=environment())
		  object$heights <- rep(1, ncol(snpset))
		  if(N > 10){
			  object$alternate.xaxis.side <- TRUE
			  side <- c(1, 3)[rep(1:2, N/2 + 1)]
			  side <- side[1:N]
			  options(warn=-1)
			  names(side) <- chromosomeNames
			  object$xaxis.side <- side
		  }
		  m <- matrix(1:(S*N), nc=N, byrow=FALSE)
		  w <- chromosomeAnnotation[paste("chr", chromosomeNames, sep=""), 1]
		  object$widths <- w/min(w)
		  object$mat <- m

                  ######################################################################
		  ##graphical parameters
		  if("copyNumber" %in% ls(assayData(snpset)) | "ratio" %in% ls(assayData(snpset))){
			  if(min(copyNumber(snpset), na.rm=TRUE) > 0) object$log <- "y" else object$log <- ""
			  ##by default, we use the same ylimit on all the plots
			  
			  ##could make plot specific by adding an option one.ylim (or
			  ##something to that effect) and calculating ylim in .plot()
			  object$ylim <- .calculateYlim(snpset, object)
		  }
		  ##---------------------------------------------------------------------------
		  ##By default, plot the cytoband at the bottom
		  ##---------------------------------------------------------------------------
		  object$cytoband.ycoords <- c(object$ylim[1], object$ylim[1] + object$cytoband.height)

		  ##---------------------------------------------------------------------------
		  ##By default, plot the predictions above the
		  ##cytoband with a little bit of space between
		  ##---------------------------------------------------------------------------
		  if(!object$add.cytoband){
			  object$hmm.ycoords <- object$cytoband.ycoords
		  } else {
			  y0 <- object$cytoband.ycoords[2] + 0.1
			  object$hmm.ycoords <- c(y0, y0 + object$height.predict)
		  }
		  
		  def.op <- options(warn=-1)
		  object$firstChromosome <- chromosomeNames[1]
		  options(def.op)

		  if(object$use.chromosome.size){
			  object$xlim <- matrix(NA, nrow=length(chromosomeNames), ncol=2)    
			  object$xlim[, 1] <- rep(0, nrow(object$xlim))
			  object$xlim[, 2] <- chromosomeSize(chromosomeNames)
			  rownames(object$xlim) <- chromosomeNames    
		  } else{
			  objList <- split(snpset, chromosome(snpset))
			  objList <- objList[chromosomeNames]
			  object$xlim <- t(sapply(objList, function(snpset) range(position(snpset))))
		  }
		  object@snpset <- snpset
		  object
		  ##set up defaults according to number of samples, chromosomes, position, etc.
	  })

setMethod("$", "ParESet", function(x, name){
	eval(substitute(snpPar(x)$NAME_ARG, list(NAME_ARG=name)))
})

setReplaceMethod("$", "ParESet",
                 function(x, name, value) {
			 snpPar(x)[[name]] = value
			 x
		 })

setMethod("snpPar", "ParESet", function(object) object@snpPar)

setReplaceMethod("snpPar", "ParESet", function(object, value) {
	object@snpPar <- value
	object
})
setMethod("snpset", "ParESet", function(object) object@snpset)
setMethod("hmmPredict", "ParESet", function(object) object@hmmPredict)

##setMethod("show", "ParESet", function(object){
##	cat("Object of class ", class(object), "\n")
##	print(snpPar(object)[1:5])
##	cat("... \n")
##})

setMethod("allPlots", "ParESet", function(object){
	list(col.axis=object$col.axis,
	     cex.main=object$cex.main,
	     cex.lab=object$cex.lab,
	     bty=object$bty,
	     ann=object$ann,
	     oma=object$oma,
	     mar=object$mar,
	     las=object$las,
	     lab=object$lab)  
})



