\name{latticeFigs}
\alias{latticeFigs}
\title{
  Generate trellis objects of log R ratios and B allele frequencies
}
\description{
  Generate trellis objects of log R ratios and B allele frequencies
}
\usage{
latticeFigs(gr, data, colors, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gr}{
    A \code{GRanges} object
}
  \item{data}{
    A \code{SummarizedExperiment} with assays "lrr" and "baf"
}
  \item{colors}{
    Colors for copy number states
}
  \item{\dots}{
    Additional arguments to \code{panel.xyplot}
}
}
\details{

}
\value{
  A list (length 2) of trellis objects with names 'lrr' and 'baf'.
}

\author{
R. Scharpf
}
\seealso{
}
\examples{
}

\keyword{dplot}

