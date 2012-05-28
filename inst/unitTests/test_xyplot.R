test_coercion <- function(){
	library(VanillaICE)
	library(oligoClasses)
	data(hmmResults, package="VanillaICE")
	data(oligoSetExample, package="oligoClasses")
	## coerce from RangedDataHMM
	dataFrame <- SNPchip:::dataFrameFromRange(hmmResults[1,], object=oligoSet)
	checkEquals(996L, nrow(dataFrame))
	## coerce from GRanges
	gr <- oligoClasses:::coerceToGRanges(hmmResults, build="hg19")
	dataFrame2 <- SNPchip:::dataFrameFromRange(gr[1,], object=oligoSet)
	checkEquals(dataFrame, dataFrame2)
}
test_xyplot <- function(){
	library(oligoClasses)
	library(IRanges)
	library(VanillaICE)
	data(oligoSetExample, package="oligoClasses")
	data(hmmResults, package="VanillaICE")
	fig <- xyplot(cn~x | range, data=oligoSet,
		      range=hmmResults,
		      frame=2e6, panel=xypanel,
		      cex=2,
		      pch=".",
		      col.het="salmon",
		      fill.het="salmon",
		      col.hom="royalblue",
		      fill.hom="royalblue",
		      state.cex=0.5,
		      border="orange", scales=list(x="free"),
		      par.strip.text=list(cex=0.5),
		      xlab="Mb", ylab=expression(log[2]("copy number")))
	checkTrue(is(fig, "trellis"))
}
