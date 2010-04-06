setClass("HmmPredict", contains="eSet",
	 representation(states="character",
			breakpoints="data.frame"))
setClassUnion("NULLorHmmPredict", c("NULL", "HmmPredict"))
setClass("ParESet",
	 representation(snpPar="list",
			snpset="SnpSet",
			hmmPredict="NULLorHmmPredict"))

setClass("ParSnpCopyNumberSet", contains="ParESet")
setClass("ParSnpCallSet", contains="ParESet")
setClass("ParSnpSet", contains="ParSnpCopyNumberSet")

setValidity("ParESet", function(object){
	valid <- validObject(object@snpset)
	if(!valid) msg <- "snpset is not a valid object" else msg <- NULL
	valid <- valid && validObject(object@hmmPredict)
	if(!valid) msg <- c(msg, "hmmPredict is not a valid object") 
	return(msg)
})

