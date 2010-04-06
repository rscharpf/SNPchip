\name{snpPar}
\alias{snpPar}
\alias{snpPar<-}
\title{Accessor for graphical parameters}
\description{
Accessor for the list of graphical parameters in objects of class \code{ParESet}.
}
\usage{
snpPar(object)
}
\arguments{
  \item{object}{object of class \code{ParESet}}
}
\value{
  \code{list}
}
\author{RS}
\seealso{ \code{\link{par}}, \code{\link[graphics]{layout}}}
\examples{
data(sample.snpset)
options(error=recover)
chr23=sample.snpset[chromosome(sample.snpset) == 23, ]
object <- plot(sample.snpset[chromosome(sample.snpset) == 23, ])
str(snpPar(object))
\dontrun{
show(object)
}
}
\keyword{methods}

