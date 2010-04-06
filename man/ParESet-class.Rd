\name{ParESet-class}
\docType{class}
\alias{ParESet-class}
\alias{$,ParESet-method}
\alias{$<-,ParESet-method}
%\alias{allPlots,ParESet-method}
%\alias{getPar,ParESet,SnpLevelSet-method}
\alias{hmmPredict,ParESet-method}
\alias{initialize,ParESet-method}
\alias{plotSnp,ParESet,SnpLevelSet-method}
\alias{snpPar,ParESet-method}
\alias{snpPar<-,ParESet-method}
\alias{snpset,ParESet-method}


\title{Class "ParESet"}
\description{A class containing a list of graphical parameters for
  plotting objects extending \code{SnpLevelSet}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ParESet", layout, col.axis, cex.main, cex.axis, cex.legend, cex, cex.lab, pch, col, bg, xaxs, xaxt, yaxs, yaxt, lab, adj, bty, ann, useLayout, mar, oma, las, log, ylab, side.ylab, outer.ylab, line.ylab, cex.ylab, xlab, outer.xlab, side.xlab, cex.xlab, line.xlab, outer.axis, line.axis, main, col.centromere, border.centromere, xlim, ylim, one.ylim, add.cytoband, outer.cytoband, outer.cytoband.axis, label.cytoband, use.chromosome.size, label.chromosome, line.label.chromosome, xaxis.side, alternate.xaxis.side, mat, heights, widths, respect, firstChromosome, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{snpPar}:}{Object of class \code{list}}
    \item{\code{hmmPredict}:}{Object of class \code{list}}
    \item{\code{snpset}:}{Object inherited from class \code{SnpLevelSet}}    
    }
}
\section{Methods}{
  \describe{
%    \item{allPlots}{\code{signature(object = "ParESet")}}
%    \item{getPar}{\code{ParESet, SnpLevelSet}: See also
%      \code{\link{getPar}}}
    \item{hmmPredict}{\code{signature(object="ParESet")}: Accessor for
     \code{HmmPredict} instance.  See also \code{\link{hmmPredict}}}
    \item{initialize}{\code{signature(.Object = "ParESet")}}
    \item{plotSnp}{\code{ParESet, SnpLevelSet}:See also \code{\link{plotSnp}}}
    \item{snpPar}{\code{signature(object = "ParESet")}}
    \item{snpPar<-}{\code{signature(object = "ParESet")}}
    \item{snpset}{\code{signature(object="ParESet")}: Accessor for
      SNP data (\code{SnpLevelSet} instance).  See \code{\link{snpset}}.}
  }
}
\author{R. Scharpf}

\seealso{
  \code{\link{ParSnpCallSet-class}},
  \code{\link{ParSnpCopyNumberSet-class}},
  \code{\link{ParSnpSet-class}}  
}
\examples{
showClass("ParESet")
}
\keyword{classes}
