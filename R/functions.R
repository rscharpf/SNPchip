centromere <- function(chromosome, build="hg18", verbose=FALSE){
	if(verbose)  message(paste("centromere coordinates based on build", build))
	if(missing(chromosome) | !(chromosome %in% c(1:22, "X"))) stop("must specify chromosome 1-22, or X as character string")
	pathto <- system.file(build, package="SNPchip")
	tmp <- read.table(file.path(pathto, "centromeres.txt"), as.is=TRUE)
##	data(chromosomeAnnotation, package="SNPchip", envir=environment())
	as.integer(tmp[paste("chr", chromosome, sep=""), ])
}

##chromosome2integer <- function(chrom){
##	chrom[chrom == "X"] <- 23; chrom[chrom == "Y"] <- 24; chrom[chrom == "XY"] <- 25; chrom[chrom=="M" | chrom == "MT"] <- 26
##	as.integer(chrom)
##}

integer2chromosome <- function(chrom){
	chrom[chrom == 23] <- "X"; chrom[chrom == 24] <- "Y"; chrom[chrom == 25] <- "XY"; chrom[chrom==26] <- "M"
	chrom
}

##chromosome2numeric <- function(chromosome){
##	chrom <- as.character(chromosome)
##	chrom[chrom == "X"] <- 23
##	chrom[chrom == "XY"] <- 24
##	chrom[chrom == "Y"] <- 25
##	chrom[chrom == "M"] <- 26
##	chrom <- as.numeric(chrom)
##	chrom
##}

chromosomeSize <- function(chromosome, build="hg18", verbose=FALSE){
	if(verbose) message(paste("chromosome size using build", build))
	if(!is.character(chromosome)) stop("argument to chromosomeSize must be one of the following character strings: 1, ..., 22, X, or Y")
	if(length(grep("chr", chromosome)) == 0) chromosome <- paste("chr", chromosome, sep="")
	if(any(!(chromosome %in% paste("chr", c(1:22, "X", "Y", "XY", "M"), sep="")))) stop("chromosome must be chr1-22, chrX, chrY, or chrM")
	pathto <- system.file(build, package="SNPchip")
	tmp <- read.table(file.path(pathto, "chromInfo.txt"), as.is=TRUE, row.names=1)
	##data(chromosomeAnnotation, package="SNPchip", envir=environment())
	tmp[chromosome, 1]
}


.labelChromosome <- function(object, op, j){
  if(j == 1){
    if(op$label.chromosome)
      mtext(unique(chromosome(object)), side=3, outer=FALSE, line=op$line.label.chromosome, cex=op$cex.lab)
  }
}

.isHomozygous <- function(object){
  calls(object) == 1 | calls(object) == 3
}

showSummary <- function(object, where, bty, legend.panel, cex, col, digits){
	f <- function(x, where, bty, legend.panel, cex, col){
		par(mar=rep(0,4))
		if(legend.panel) plot(0:1, 0:1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
		legend(where, legend=c(
			      paste(x["ho"], " %AA/BB", sep=""),
			      paste(x["ht"], " %AB", sep=""),
			      paste(x["cn"], " avg CN"),
			      paste(x["cn.sd"], " sd")), bty=bty,
		       title=substr(x["samplenames"], 1, 10),
		       y.intersp=1.1,
		       cex=cex,
		       text.col=c("black", col[1], col[2], "black", "black"))
	}
	if(sum(chromosome(object) != "X") > 0){
		obj <- object[chromosome(object) != "X", ]
	} else { obj <- object}
	cn <- colMeans(as.matrix(copyNumber(obj)))
	cn.sd <- apply(copyNumber(obj), 2, sd)
	ht <- colMeans(ifelse(calls(obj) == 2, 1, 0))
	ho <- colMeans(ifelse(calls(obj) == 1 | calls(obj) == 3, 1, 0))
	stats <- cbind(cn, cn.sd, ht, ho)
	stats <- round(stats, digits)
	stats <- data.frame(stats); stats$samplenames <- sampleNames(object)
	apply(stats, 1, f, where=where, bty=bty,
	      legend.panel=legend.panel, cex=cex, col=col)
}







