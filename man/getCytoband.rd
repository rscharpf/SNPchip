\name{getCytoband}
\alias{getCytoband}
\title{getCytoband}

\description{
This function generates a \code{data.frame} with the respective cytoband names, chromosomes, Giemsa stain, and the start and end positions.   These tables can then be used to plot chromosome idiograms.  Currently, cytoband annotation for UCSC genome builds hg18 and hg19 are supported.

\usage{
  getCytoband(build)

\arguments{
\item{build} {A character object or vector (e.g, "hg19").}
}


\value{
\code{data.frame}
}

\examples{
cytoband <- getCytoband("hg19")
cytoband <- cytoband[cytoband$chr == "chr1", ]
plotIdiogram(1, "hg18", cytoband=cytoband, cex.axis=0.6)
}

\author{Michael "Destroyer of Worlds" Considine, Rob "The Boy Wonder" Scharpf}



\seealso{\code{\link{plotIdiogram}}}
\keyword{misc}