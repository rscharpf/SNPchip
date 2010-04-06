\name{featureData-accessors}
\alias{chromosome}
\alias{chromosome-methods}
\alias{chromosome,SnpLevelSet-method}
\alias{chromosome,SnpCallSet-method}
\alias{chromosome,SnpCopyNumberSet-method}
\alias{chromosome,oligoSnpSet-method}
\alias{chromosome,AnnotatedDataFrame-method}

\alias{position}
\alias{position-methods}
\alias{position,SnpLevelSet-method}
\alias{position,SnpCallSet-method}
\alias{position,SnpCopyNumberSet-method}
\alias{position,oligoSnpSet-method}
\alias{position,AnnotatedDataFrame-method}
\alias{position,SnpLevelSet-method}

\title{Accessors for feature-level SNP annotation}

\description{
  
  \code{featureData} accessors for classes defined in SNPchip

}


\arguments{

  \item{object}{object inheriting from \code{SnpLevelSet}}

}

\details{
  \code{chromosome} extracts character string of chromosome number for
  each SNP in the \code{object}.

  \code{position} extracts the physical position (base pair number) for
  each SNP in the \code{object}.

}

\author{R. Scharpf}
\keyword{methods}
\keyword{manip}

