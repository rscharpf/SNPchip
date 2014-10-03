\name{dataFrame}
\alias{dataFrame}

\title{
  Generic function for coercing gSet objects to data.frame
}

\description{
  Generic function for coercing \code{gSet} objects to \code{data.frame} as a
  precursor to plotting with lattice
}

\usage{
dataFrame(range, data, ...)
}

\arguments{

  \item{range}{
    A \code{GenomicRanges} object containing interval(s) for which
  low-level data should be plotted
}

  \item{data}{ A container for the low-level data (e.g.,
    \code{BafLrrSet}) or a \code{SummarizedExperiment} }

  \item{\dots}{
    Additional arguments passed to \code{findOverlaps}. E.g., argument
    \code{maxgap} can be used to select the size of the window
    surrounding the genomic intervals in \code{range} for plotting.
  }

}

\value{
  A \code{data.frame} with column labels that depend on the class of \code{data}.
}

\author{
R. Scharpf
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{methods}% __ONLY ONE__ keyword per line
