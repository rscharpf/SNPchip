\name{ParSnpCopyNumberSet-class}
\docType{class}
\alias{ParSnpCopyNumberSet-class}
\alias{initialize,ParSnpCopyNumberSet-method}

\title{Class "ParSnpCopyNumberSet" }
\description{Graphical parameters for plotting \code{SnpCopyNumberSet} objects}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParSnpCopyNumberSet", ...)}.
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
    \item{initialize}{\code{signature(.Object = "ParSnpCopyNumberSet")}}
    \item{plotSnp}{\code{signature(object = "ParSnpCopyNumberSet", snpset = "SnpCopyNumberSet")}}
  }
}
\author{R. Scharpf}
\examples{
showClass("ParSnpCopyNumberSet")
}
\keyword{classes}
