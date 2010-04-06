\name{.calculateYlim}
\alias{.calculateYlim}
\title{Computes y-limits of plotting device}
\description{
  Computes y-limits of plotting device
}
\usage{
.calculateYlim(object, op)
}
\arguments{
  \item{object}{object inheriting from \code{eSet}}
  \item{op}{object inherting from \code{ParESet}}
}
\details{
  
  If element \code{copyNumber} is present in the \code{assayData}, the
  range for the y-axis on the plotting device is calculated as the .001
  and 0.999 quantiles of the copyNumber.

  If element \code{copyNumber} is not present, we plot the genotype
  calls on the vertical axis.
  
}
\value{
  \code{numeric}
}
\references{}
\author{}
\note{
}
\seealso{}
\examples{
}
\keyword{manip}
