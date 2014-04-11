\name{centromere}
\alias{centromere}
\title{Coordinates of centromere}
\description{
  Extracts coordinates of centromere for a particular chromosome
}
\usage{
centromere(chromosome, build, verbose=FALSE)
}
\arguments{
  \item{chromosome}{Chromosome name. Several formats for specifying
    chromosome are allowed (see examples). }
  \item{build}{character string. Supported UCSC builds are `hg18' and `hg19'.}
  \item{verbose}{Logical.  Displays build used to annotate the
    centromere coordinates when TRUE}
}
\value{
  integer: start and stop coordinates of centromere in basepairs
}
\author{R. Scharpf}
\examples{
centromere(1, "hg18")
centromere("1", "hg18")
centromere("chr1", "hg18")
centromere(1, "hg19")
centromere("X", "hg18")
}
\keyword{manip}

