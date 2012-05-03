\name{xyplotLrrBaf}
\alias{xyplotLrrBaf}
\title{xyplot lattice function for RangedData and oligoSnpSet objects}
\description{
  For each genomic interval in the ranged data, a plot of the log R ratios and B allele frequencies stored in
  the oligoSnpSet are plotted.
}
\usage{
xyplotLrrBaf(rd, object, frame, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rd}{
    A \code{RangedDataCNV}-derived object.
}
  \item{object}{
    A \code{oligoSnpSet} or \code{BeadStudioSet} object with assayData
    elements for log R ratios and B allele frequencies.
  }

  \item{frame}{ The genomic distance in basepairs to plot on either side
    of the genomic interval in the \code{rd} object.  }

  \item{\dots}{Additional arguments passed to the panel function. See details.
}

}

\details{ The \code{xypanelBaf} function is a panel function that does
  the actual plotting of the genomic data.

  }

\value{
  A \code{trellis} object.
}

\author{
R. Scharpf
}

\seealso{
  \code{\link{xypanelBaf}}
}
\examples{
library(crlmm)
library(VanillaICE)
data(cnSetExample, package="crlmm")
oligoSetList <- cnSet2oligoSetList(cnSetExample)
fit <- hmm(oligoSetList, p.hom=0)
## We're interested in this range
rangeOfInterest <- RangedDataCNV(IRanges(3.7e6, 5.9e6), chrom=8, sampleId="NA19007")
index <- subjectHits(findOverlaps(rangeOfInterest, fit))
xyplotLrrBaf(fit[index, ], oligoSetList[[1]], frame=1e6, panel=xypanelBaf, cex=0.2)
}
\keyword{hplot}

