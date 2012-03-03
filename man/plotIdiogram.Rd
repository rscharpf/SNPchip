\name{plotIdiogram}
\alias{plotIdiogram}
\title{Plots idiogram for one chromosome}

\description{ Draw an idiogram for the specified chromosome.}

\usage{
plotIdiogram(chromosome, cytoband, cytoband.ycoords, xlim, ylim=c(0, 2),
new=TRUE, label.cytoband=TRUE, label.y=NULL, srt, cex.axis=1,
outer=FALSE, taper=0.15, verbose=FALSE, build="hg18", ...)
}
\arguments{
  \item{chromosome}{character string or integer: which chromosome to draw the
    cytoband}
  \item{cytoband}{data.frame containing cytoband information}
  \item{cytoband.ycoords}{numeric: y coordinates }
  \item{xlim}{x-axis limits}
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
  \item{build}{Character string.  Currently only "hg18" is allowed.}
  \item{\dots}{additional arguments to plot}
}

\author{Robert Scharpf and Jason Ting}
\examples{
plotIdiogram("1")
}
\keyword{aplot}

