\name{hapmap}
\alias{hapmap}
\docType{data}
\title{Data from the hapmap project}
\description{
Data obtained from the hapmap project and pre-processed by CRLMM to
produce  genotype calls.  Due to size
limitations, we only included every 10th SNP from the Affy's 50k Xba
chip.  'Copy number' estimates were obtained from using  the calculate
copy number function in the oligoClasses package and is used solely to
illustrate properties of the classes defined in SNPchip.  Better approaches
for estimating copy number are available in other packages, and
currently under development in the oligo package.
}
\usage{data(hapmap)}
\author{R. Scharpf}
\examples{
data(hapmap)
}
\keyword{datasets}
