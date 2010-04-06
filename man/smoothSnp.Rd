\name{smoothSnp}
\alias{smoothSnp}
\alias{smoothSnp,oligoSnpSet-method}
\title{A simple nonparametric smoother for genotype and copy number}
\description{
  A loess smoother for finding regions of reduced copy number and loss
  of heterozygosity.
}
\usage{
smoothSnp(object, ...)
}
\arguments{
  \item{object}{Object derived from \code{AnnotatedSnpSet}}
  \item{\dots}{For additional arguments see details.}
}
\details{
  
  Additional arguments include the chromosomes (character vector) and
  samples (numerical vector) to smooth.  Options for smoothing are set
  by specifying the span and the method.  See the SNPchip vignette.

  See the R package VanillaICE for a hidden Markov model for more formal
  inference regarding regions of LOH and copy number alterations.
  
}

\value{

  An object of the same class, e.g., \code{AnnotatedSnpSet}, where
  assayData elements copyNumber and calls are replaced by the smoothed
  values.

}

\references{RS}
\author{Robert Scharpf}
\keyword{methods}

