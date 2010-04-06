\name{centromere}
\alias{centromere}
\title{Coordinates of centromere}
\description{
  Extracts coordinates of centromere for a particular chromosome
}
\usage{
centromere(chromosome, build="hg18", verbose=FALSE)
}
\arguments{
  \item{chromosome}{character string: "1", ..., "22", "X", or "Y"}
  \item{build}{character string. Currently only build 'hg18' is
    provided}
  \item{verbose}{Logical.  Displays build used to annotate the
    centromere coordinates when TRUE}
}
\value{
  integer: start and stop coordinates of centromere in basepairs
}
\author{R. Scharpf}
\examples{
centromere("1")
}
\keyword{manip}

