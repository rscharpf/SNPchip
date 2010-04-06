\name{plotSnp}
\alias{plotSnp}
\alias{plot}
\alias{plot,eSet-method}
\alias{plotSnp,eSet-method}
\alias{plotSnp,ParSnpCopyNumberSet,SnpCopyNumberSet-method}
\alias{plotSnp,ParSnpCallSet,SnpCallSet-method}
\alias{plotSnp,ParSnpSet,oligoSnpSet-method}
\alias{plotSnp,SnpLevelSet-method}
\alias{plot,ParSnpCopyNumberSet,SnpCopyNumberSet-method}
\alias{plot,ParSnpCallSet,SnpCallSet-method}
\alias{plot,ParSnpSet,oligoSnpSet-method}
\alias{plot,SnpLevelSet-method}
\title{Plots copy number and genotype calls against physical position }
\description{

  Returns an object inheriting from class \code{ParESet} -- essentially,
  a list of default graphical parameters that can be modified as needed.
  The \code{show} method for this class plots the copy number and
  genotype calls versus physical position for an arbitrary number of
  samples and chromosomes.

}
\usage{
plotSnp(object, hmmPredict, ...)
\S4method{plot}{SnpLevelSet}(x, y, ...)
}

\arguments{
  \item{object, x}{An object extending \code{ParESet}}
  \item{hmmPredict, y}{An object of class \code{HmmPredict}}
  \item{\dots}{ Additional arguments to the initialization methods of
  the \code{ParESet} classes}
}
\details{
  See examples in the vignette
}

\author{R. Scharpf}

\seealso{
  \code{\link{par}},
  \code{\link{ParESet-class}},
  \code{\link{ParSnpCallSet-class}},
  \code{\link{ParSnpCopyNumberSet-class}},
  \code{\link{ParSnpSet-class}}
}

\examples{
data(sample.snpset)
graph.par <- plot(sample.snpset[chromosome(sample.snpset) < 5, ])

class(graph.par)
\dontrun{
	## to plot the data, the following are equivalent
	graph.par
	# or
	show(graph.par)
	# or
	print(graph.par)
}
}
\keyword{hplot}

