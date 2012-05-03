\name{xypanelBaf}
\alias{xypanelBaf}

\title{Panel function for plotting copy number and B allele frequencies
  for a genomic interval.}
\description{Panel function for plotting copy number and B allele frequencies
  for a genomic interval.
}
\usage{
xypanelBaf(x, y, gt, baf, is.snp, range, col.hom = "grey20", fill.hom = "lightblue", col.het = "grey20", fill.het = "salmon", col.np = "grey20", fill.np = "grey60", show.state = TRUE, state.cex = 1, col.state = "blue", ..., subscripts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    physical position in basepairs
  }

  \item{y}{
    total copy number (relative or absolute)
  }

  \item{gt}{
    Genotypes coded as integers (1=AA, 2=AB, 3=BB).  This is optional.
    If provided one can color code the plotting symbols by the genotype.
  }

  \item{baf}{
    B allele frequencies.
}

  \item{is.snp}{
    Logical.  Indicator of whether the marker hybridized to a known SNP
    or a nonpolymorphic region of the genome.
  }

  \item{range}{
    A \code{RangedDataCNV}-derived object indicating the genomic interval to plot.
}

  \item{col.hom}{
    Color to use for homozygous genotypes.
  }

  \item{fill.hom}{
    Fill color to use for homozygous genotypes.
  }

  \item{col.het}{
    Color to use for heterozygous genotypes.
}
  \item{fill.het}{
    Fill color to use for heterozygous genotypes.
}
  \item{col.np}{
    Color to use for nonpolymorphic markers
}
  \item{fill.np}{
    Fill color for nonpolymorphic markers.
}
  \item{show.state}{
    Logical indicating whether to display the copy number state for a \code{RangedDataHMM} object.
}
  \item{state.cex}{
    Size of the font for displaying the HMM state. Ignored if
    \code{show.state} is FALSE.
}
  \item{col.state}{
    Color for displaying the state.
}
  \item{\dots}{
    Additional arguments passed to \code{panel.xyplot}.
}
  \item{subscripts}{
    See \code{panel.xyplot}
}
}
\details{
  Function for plotting B allele frequencing and copy number on a
    trellis display.  Intended to be passed to the panel argument of the
    function \code{xyplotLrrBaf} and should not be called directly by
    the user.
}


\author{
R.Scharpf
}
\seealso{
  \code{\link{xyplotLrrBaf}}
}
\keyword{aplot}

