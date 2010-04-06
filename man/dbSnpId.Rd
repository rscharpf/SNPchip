\name{dbSnpId}
\alias{dbSnpId}
\alias{dbSnpId,AnnotatedDataFrame-method}
\alias{dbSnpId,SnpLevelSet-method}
\title{The dbSNP identifier }
\description{
  Unique identifier for each SNP assigned by NIH
}
\usage{
dbSnpId(object)
}
\arguments{
  \item{object}{Object possibly derived from \code{eSet} or
    \code{AnnotatedDataFrame} class.}
}
\value{
  A (typically very long) character vector of SNP identifiers.
}
\author{R. Scharpf}

\examples{
\dontrun{
data(sample.snpset)
dbSnpId(sample.snpset)[1:5]
}
}
\keyword{manip}

