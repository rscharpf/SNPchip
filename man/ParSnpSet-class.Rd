\name{ParSnpSet-class}
\docType{class}
\alias{ParSnpSet-class}
\alias{initialize,ParSnpSet-method}

\title{Class "ParSnpSet"}
\description{Graphical parameters for \code{oligoSnpSet} objects}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParSnpSet", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{snpPar}:}{Object of class \code{"list"}}
  }
}
\section{Extends}{
Class \code{\linkS4class{ParSnpCopyNumberSet}}, directly.
Class \code{\linkS4class{ParESet}}, by class "ParSnpCopyNumberSet", distance 2.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "ParSnpSet")}}
    \item{plotSnp}{\code{signature(object = "ParSnpSet", snpset = "oligoSnpSet")}}
  }
}
\author{R. Scharpf}

\seealso{
  \code{\link[oligoClasses:oligoSnpSet-methods]{oligoSnpSet-class}}
}
\examples{
showClass("ParSnpSet")
}
\keyword{classes}
