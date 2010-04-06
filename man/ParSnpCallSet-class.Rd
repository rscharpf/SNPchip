\name{ParSnpCallSet-class}
\docType{class}
\alias{ParSnpCallSet-class}
\alias{initialize,ParSnpCallSet-method}

\title{Class "ParSnpCallSet"}
\description{Graphical parameters for object of class \code{SnpCallSet}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParSnpCallSet", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{snpPar}:}{Object of class \code{"list"} }
  }
}
\section{Extends}{
  Class \code{\linkS4class{ParESet}}, directly.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "ParSnpCallSet")} }
    \item{plotSnp}{\code{signature(object = "ParSnpCallSet", snpset = "SnpSet")}}
  }
}
\author{R.Scharpf}
\examples{
showClass("ParSnpCallSet")
}
\keyword{classes}
