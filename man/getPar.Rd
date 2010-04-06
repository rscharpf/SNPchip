\name{getPar}
\alias{getPar}
\alias{getPar,ParESet-method}
\title{Adds graphical parameters for plotting SNP data}
\description{
  Not intended to be called directly.
  Adds graphical parameters for plotting SNP data to one of the
  graphical parameter classes (e.g., \code{ParSnpSet},
  \code{ParSnpCopyNumberSet}, etc).  
}

\usage{
  getPar(object, ...)
}

\arguments{
  \item{object}{ an object inheriting from \code{ParESet}}
  \item{\dots}{additional arguments to \code{par}}
}
\details{

Adds graphical parameters to an object inherited from class
      \code{ParESet} that depend on the data class object (an object
      inherited from class \code{SnpLevelSet}). For instance, graphical
      parameters specifying the layout depend on the number of samples
      and chromosomes in the object \code{snpset}
      
}

\value{
  An object of the same class as \code{object}
}

\author{R. Scharpf}
\seealso{plotSnp}
\keyword{manip}
\keyword{methods}
