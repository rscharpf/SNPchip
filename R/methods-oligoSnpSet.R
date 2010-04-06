setMethod("show", "oligoSnpSet", function(object) {
	cat(class( object ), " (storageMode: ", storageMode(object), ")\n", sep="")
	adim <- dim(object)
	if (length(adim)>1)
		cat("assayData:",
		    if (length(adim)>1)
		    paste(adim[[1]], "features,",
			  adim[[2]], "samples") else NULL,
		    "\n")
	cat("  element names:",
	    paste(assayDataElementNames(object), collapse=", "), "\n")
	cat("experimentData: use 'experimentData(object)'\n")
	pmids <- pubMedIds(object)
	if (length(pmids) > 0 && all(pmids != ""))
		cat("  pubMedIds:", paste(pmids, sep=", "), "\n")
	cat("Annotation:", annotation(object), "\n")
	cat("phenoData\n")
	if(length(adim) > 1) show(phenoData(object))
	cat("featureData\n")      
	if(length(adim) > 1) show(featureData(object))
	cat("Annotation ")
  show(annotation(object))
})


setMethod("summary", "oligoSnpSet",
          function(object, digits=3, MISSING.CODE=4, ...){
            ##Calculate mean,sd copy number and prop no calls, prop het calls
            ##for each chromosome in each sample.  Return an S x 23 matrix for
            ##each.
		  require(genefilter) || stop("package genefilter is not available")
		  obj <- split(object, chromosome(object))
		  means <- do.call(rbind, lapply(obj, function(x) colMeans(copyNumber(x))))

		  colSds <- function(x) rowSds(t(x))            
		  sds <- do.call(rbind, lapply(obj, function(x) colSds(copyNumber(x))))
            
		  ##Proportion of no calls for each chromosome
		  if(!is.na(MISSING.CODE)){
			  pNoCall <- do.call(rbind, lapply(obj, function(x, MISSING.CODE) colMeans(calls(x) == MISSING.CODE), MISSING.CODE))
		  } else{
			  pNoCall <- do.call(rbind, lapply(obj, function(x, MISSING.CODE) colMeans(is.na(calls(x))), MISSING.CODE))
		  }
		  pHom <- do.call(rbind, lapply(obj, function(x) colMeans(calls(x) == 1 | calls(x) == 3)))
		  pHet <- do.call(rbind, lapply(obj, function(x) colMeans(calls(x) == 2)))
		  x <- list(means, sds, pHom, pHet, pNoCall)
		  names(x) <- c("avg.CN", "sd.CN", "prop.Hom", "prop.Het", "prop.NoCall")            
		  x <- lapply(x, round, digits)
		  return(x)
          })


#create object of class snpscan with smoothed copynumbers and smoothed loh calls
#LOH calls are smoothed by setting homozygous = 0 and heterozygous = 1
setMethod("smoothSnp", "oligoSnpSet",
          function(object, 
                   span=1/10,
                   method="loess",
                   imputeNoCalls=TRUE,
                   verbose=TRUE){
            ##convert homozygous to 0 and heterozygous to 1
            calls(object)[calls(object) == 1 | calls(object) == 3] <- 0
            calls(object)[calls(object) == 2] <- 1
            smoothChromosome <- function(obj, span){
              loessX <- function(X, location, span){
                fit <- loess(X ~ location, span = span)$fitted
                return(fit)
              }
              ##Order by physical position before smoothing
              obj <- obj[order(position(obj)), ]
              cn.smooth <- apply(copyNumber(obj), 2, loessX, position(obj), span=span)
              call.smooth <- apply(calls(obj), 2, loessX, location=position(obj), span=span)
              rownames(cn.smooth) <- rownames(call.smooth) <- featureNames(obj) 
              copyNumber(obj) <- cn.smooth
              calls(obj) <- call.smooth
              obj
            }
            object.list <- split(object, chromosome(object))
            smooth.list <- lapply(object.list, smoothChromosome, span=span)
            smooth.set <- unsplitSnpSet(smooth.list,
                                        featureData(object),
                                        copyNumber=do.call(rbind, lapply(smooth.list, copyNumber)),
                                        cnConfidence=do.call(rbind, lapply(smooth.list, cnConfidence)),
                                        calls=do.call(rbind, lapply(smooth.list, calls)),
                                        callProbability=do.call(rbind, lapply(smooth.list, calls)),
                                        phenoData=phenoData(smooth.list[[1]]),
                                        annotation=annotation(smooth.list[[1]]),
                                        experimentData=experimentData(smooth.list[[1]]))
            return(smooth.set)
          })


## setAs("oligoSnpSet", "SnpSet",
##       function(from, to){
## 	      new("SnpSet",
## 		  call=calls(from),
## 		  callProbability=assayData(from)[["callProbability"]],
## 		  featureData=featureData(from),
## 		  experimentData=experimentData(from),
## 		  phenoData=phenoData(from),
## 		  annotation=annotation(from),
## 		  protocolData=protocolData(from))
##       })

