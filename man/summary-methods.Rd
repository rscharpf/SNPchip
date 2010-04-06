\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,SnpSet-method}
\alias{summary,AnnotatedSnpSet-method}
\alias{summary,oligoSnpSet-method}
\alias{summary,eSet-method}
\title{Summary statistics for various SNP classes that extend eSet}
\description{
Methods for function \code{summary} 
}
\section{Methods}{
\describe{

\item{object = "SnpSet"}{calculates percent homozygous calls, and
  percent heterozygous calls for each chromosome, and then aggregated
  over all autosomes.}

\item{object = "AnnotatedSnpSet"}{calculates average copy number,
  standard deviation of copy number, percent homozygous calls, and
  percent heterozygous calls for each chromosome, and then aggregated
  over all autosomes.}
  }} 
\keyword{methods}

