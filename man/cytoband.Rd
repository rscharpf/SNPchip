\name{cytoband}
\alias{cytoband}
\docType{data}
\title{Start and stop sites of cytoband}
\description{
Contains start and stop sites of cytoband for Build 35
}
\usage{data(cytoband)}
\format{
  A data frame with 862 observations on the following 5 variables.
  \describe{
    \item{\code{chrom}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr13} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr21} \code{chr22} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrX} \code{chrY}}
    \item{\code{chromStart}}{a numeric vector}
    \item{\code{chromEnd}}{a numeric vector}
    \item{\code{name}}{a factor indicating which chromosomal arm}
    \item{\code{gieStain}}{a factor with levels \code{acen} \code{gneg} \code{gpos100} \code{gpos25} \code{gpos50} \code{gpos75} \code{gvar} \code{stalk}}
  }
}

\source{
\url{http://pevsnerlab.kennedykrieger.org/snpscan_07_sourcecode.htm}
}

\examples{
data(cytoband)
}
\keyword{datasets}
