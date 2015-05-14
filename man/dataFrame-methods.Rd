\name{dataFrame-methods}
\docType{methods}
\alias{dataFrame-methods}
\alias{dataFrame,GRanges,gSet-method}
\alias{dataFrame,GRanges,RangedSummarizedExperiment-method}
\title{Construct a data.frame from genomic data for plotting}
\description{
  Construct a \code{data.frame} of genomic data (log R ratios and BAFs)
  from a \code{RangedSummarizedExperiment} with markers in the interval given
  by the \code{GRanges} object.
}
\section{Methods}{
\describe{

\item{\code{signature(range = "GRanges", data = "gSet")}}{ The argument
  \code{range} is often intervals from a hidden Markov model fit to the
  genomic data in the \code{data} object. \code{gSet}-derived classes
  contain assay data on copy number and allele frequencies.
}

\item{\code{signature(range = "GRanges", data =
  "RangedSummarizedExperiment")}}{ The argument \code{range} is often
  intervals from a hidden Markov model fit to the genomic data in the
  \code{data} object. The \code{RangedSummarizedExperiment} is assumed to
  contain log R ratio (lrr) and B allele frequency (baf) assays.
}
}

}

\keyword{methods}

