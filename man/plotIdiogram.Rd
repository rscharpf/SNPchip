\name{plotIdiogram}
\alias{plotIdiogram}
\alias{plotCytoband2}
\title{Plots idiogram for one chromosome}

\description{ Draw an idiogram for the specified chromosome.}

\usage{
plotIdiogram(chromosome, build, cytoband, cytoband.ycoords, xlim, ylim=c(0, 2),
new=TRUE, label.cytoband=TRUE, label.y=NULL, srt, cex.axis=1,
outer=FALSE, taper=0.15, verbose=FALSE, unit=c("bp", "Mb"),
is.lattice=FALSE,...)
plotCytoband2(chromosome, build, cytoband, xlim, xaxs="r", new=TRUE,
label.cytoband=TRUE, cex.axis=1, outer=FALSE, verbose=TRUE, ...)
}
\arguments{
  \item{chromosome}{character string or integer: which chromosome to draw the
    cytoband}
  \item{build}{UCSC genome build. Supported builds are "hg18" and "hg19".}
  \item{cytoband}{data.frame containing cytoband information}
  \item{cytoband.ycoords}{numeric: y coordinates }
  \item{xlim}{x-axis limits}
  \item{xaxs}{numeric. See \code{par}}
  \item{ylim}{y-axis limits}
  \item{new}{logical: new plotting device}
  \item{label.cytoband}{logical: if TRUE, labels the cytobands}
  \item{label.y}{numeric: height (y-coordinate) for cytoband label}
  \item{srt}{string rotation for cytoband labels.  See \code{par}}
  \item{cex.axis}{size of cytoband labels. See par}
  \item{outer}{logical: whether to draw the labels in the outer
    margins. See par}
  \item{taper}{tapering for the ends of the cytoband}
  \item{verbose}{Logical.  If TRUE, displays human genome build used to
    annotated the cytoband coordinates.}
  \item{unit}{Character string indicating the unit for physical position
    on the x-axis. Available options are basepairs (bp) or Mb.}
  \item{is.lattice}{logical indicating whether your drawing the cytoband
    on a lattice graphic.}
  \item{\dots}{additional arguments to plot}
}

\author{Robert Scharpf and Jason Ting}
\examples{
plotIdiogram("1", "hg18")
plotIdiogram("1", "hg19")
plotIdiogram("1", build="hg19", cex=0.8, label.cytoband=FALSE)
## user-defined coordinates
plotIdiogram("1", build="hg19", cex=0.8, label.cytoband=FALSE,
ylim=c(0,1), cytoband.ycoords=c(0.1, 0.3))

library(oligoClasses)
sl <- getSequenceLengths("hg19")[c(paste("chr", 1:22, sep=""), "chrX", "chrY")]
ybottom <- seq(0, 1, length.out=length(sl)) - 0.01
ytop <- seq(0, 1, length.out=length(sl)) + 0.01
for(i in seq_along(sl)){
	chr <- names(sl)[i]
	if(i == 1){
		plotIdiogram("1", build="hg19", cex=0.8, label.cytoband=FALSE, ylim=c(-0.05,1.05), cytoband.ycoords=c(ybottom[1], ytop[1]),
			     xlim=c(0, max(sl)))
	}
	if(i > 1){
		plotIdiogram(names(sl)[i], build="hg19", cex=0.8, label.cytoband=FALSE, cytoband.ycoords=c(ybottom[i], ytop[i]), new=FALSE)
	}
}
axis(1, at=pretty(c(0, max(sl)), n=10), labels=pretty(c(0, max(sl)), n=10)/1e6, cex.axis=0.8)
mtext("position (Mb)", 1, line=2)
par(las=1)
axis(2, at=ybottom+0.01, names(sl), cex.axis=0.6)
}
\keyword{aplot}

