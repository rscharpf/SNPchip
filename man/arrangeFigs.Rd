\name{arrangeFigs}
\alias{arrangeFigs}

\title{
  Arranging two trellis objects on a grid.
}
\description{
  Helper function for arranging multipanel displays of log R ratios and
  B allele frequencing in a single figure
}
\usage{
arrangeFigs(lattice.figs, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lattice.figs}{
    A named list ('lrr' and 'baf') of two \code{trellis} object.
}
  \item{\dots}{
    ignored
}
}


\value{
  nothing
}
\author{
R. Scharpf
}

\seealso{
  \code{latticeFigs}, \code{\link{arrangeSideBySide}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
