centromereData <- function(genome){
	path <- system.file("extdata", package="SNPchip")
	path2 <- system.file("extdata", package="oligoClasses")
	## hg19:  http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz
	## hg18:  http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/gap.txt.gz
	dat <- read.delim(file.path(path, paste("gap_", genome, ".txt", sep="")), header=FALSE)
	dat <- dat[dat[, 8]=="centromere", ]
	load(file.path(path2, paste("seqlengths_", genome, ".rda", sep="")))
	gr <- GRanges(as.character(dat[, 2]), IRanges(dat[,3], dat[,4]), seqlengths=seqlengths[as.character(dat[,2])])
	gr <- sort(gr)
	saveRDS(gr, file=file.path("~/Software/SNPchip/inst/extdata", paste("gr_", genome, ".rda", sep="")))
}
