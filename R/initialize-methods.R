setMethod("initialize", "ParSnpCallSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$col <- c("lightblue", "red", "lightblue")
            .Object$ylab <- "genotype call"
            .Object
          })

setMethod("initialize", "ParSnpCopyNumberSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$ylab <- "copy number"
            .Object
          })

setMethod("initialize", "ParSnpSet",
          function(.Object, ...){
            .Object <- callNextMethod(.Object, ...)
            .Object$col <- c("lightblue", "red", "lightblue")
            .Object$bg <- c("lightblue", "red", "lightblue")
            .Object$ylab <- "copy number"
	    .Object$legend.fill <- .Object$col
            .Object
          })


setMethod("initialize", "ParESet",
          function(.Object,
		   snpset,
		   hmmPredict=NULL, ...){
		  if("snpset" %in% names(list(...))){
			  snpset <- list(...)[["snpset"]]
		  }
##		  if(missing(snpset)){
##			  stop("snpset missing")
##		  }
		  if(!extends(class(snpset), "eSet")){
			  stop("snpset must extend eSet")
		  }
		  .Object@snpset <- snpset
		  if("hmmPredict" %in% names(list(...))){
			  hmmPredict <- list(...)[["hmmPredict"]]
			  if(!extends(class(hmmPredict, "eSet"))){
				  stop("hmmPredict must extend eSet. See VanillaICE package")
			  }
			  .Object@hmmPredict <- hmmPredict
		  }
##		  if(!is.null(hmmPredict)){
##			  if(!extends(class(hmmPredict, "eSet"))){
##				  stop("hmmPredict must extend eSet. See VanillaICE package")
##			  }
##		  }
		  .Object@snpPar <- list(layout=TRUE,
					 col.axis="brown",
					 cex.main=1,
					 cex.axis=1,
					 cex.legend=1,
					 cex=2,
					 cex.lab=1,
					 pch=".",
					 col="black",
					 bg="white",
					 xaxs="r",
					 xaxt="s",
					 yaxs="r",
					 yaxt="s",
					 at=1:4,
					 lab=c(2, 5, 7), ##see par
					 adj=0,
					 bty="n",
					 ann=FALSE,
					 useLayout=TRUE,
					 mar=c(0.5, 0, 0.5, 0.2),
					 oma=c(4, 4, 4, 0.5),
					 las=1,
					 log="",
					 ylab="",
					 side.ylab=2,
					 outer.ylab=TRUE,
					 line.ylab=3,
					 cex.ylab=1,
					 xlab="Mb",
					 outer.xlab=TRUE,
					 side.xlab=1,
					 cex.xlab=1,
					 line.xlab=1,
					 outer.axis=TRUE,
					 line.axis=0,
					 main="",
					 add.centromere=FALSE,
					 col.centromere="bisque",
					 border.centromere="bisque",
					 xlim=NULL,
					 ylim=NULL,
					 one.ylim=TRUE,
					 add.cytoband=TRUE,
					 outer.cytoband=FALSE,
					 outer.cytoband.axis=FALSE,
					 label.cytoband=FALSE,
					 cytoband.ycoords=NULL,
					 hmm.ycoords=NULL,
					 cytoband.srt=90,
					 cytoband.label.y=NULL,
					 cytoband.taper=0.15,
					 cytoband.height=0.1,
					 use.chromosome.size=FALSE, #for x-axis limits
					 label.chromosome=TRUE,
					 line.label.chromosome=2,
					 xaxis.side=1,
					 alternate.xaxis.side=FALSE,
					 mat=new("matrix", 1, 1),
					 heights=1,
					 widths=1,
					 respect=FALSE,
					 firstChromosome="1",
					 abline=FALSE,
					 abline.h=2,
					 abline.col="grey80",
					 abline.lty=1,
					 abline.lwd=1,
					 abline.v=FALSE,
					 abline.v.pos=NULL,
					 abline.v.col="grey20",
					 abline.v.lty=2,
					 abline.v.lwd=0.8,
					 legend=NULL,
					 legend.location="topright",
					 legend.bty="n",
					 legend.fill="white",
					 legend.pch=21,
					 border.prediction=grey(0.4),
					 legend.predict=NULL,
					 legend.location.predict="topleft",
					 legend.fill.predict=NULL,
					 col.predict=NULL,
					 height.predict=0.2,
					 cytoband.side=1)
            .Object
          })
