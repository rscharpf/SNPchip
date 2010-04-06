\name{chromosomeSize}
\alias{chromosomeSize}
\title{Size of chromosome in number of base pairs.}
\description{
  Size of chromosome in base pairs.
}
\usage{
chromosomeSize(chromosome, build="hg18", verbose=FALSE)
}

\arguments{
  \item{chromosome}{character string for the chromosome, e.g., "1", ...,
    "22", "X", or "Y"}
  \item{build}{character string: genome build.  Currently only 'hg18' is
    provided.}
  \item{verbose}{Logical.  If TRUE, build is displayed}
}

\value{
  Numerical vector.
}

\author{Robert Scharpf}
\examples{
chromosomeSize("1", verbose=TRUE)
}
\keyword{manip}

