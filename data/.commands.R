##library(SNPchip) ##old version
##data(annSnpset)
##
##sample.snpset <- new("oligoSnpSet",
##                     calls=calls(annSnpset),
##                     callsConfidence=callsConfidence(annSnpset),
##                     copyNumber=copyNumber(annSnpset),
##                     cnConfidence=cnConfidence(annSnpset),
##                     annotation=annotation(annSnpset),
##                     featureData=featureData(annSnpset),
##                     phenoData=phenoData(annSnpset),
##                     experimentData=experimentData(annSnpset),
##                     annotation="pd.mapping50k.xba240")
##save(sample.snpset, file="~/projects/software/SNPchip/data/sample.snpset.RData", compress=TRUE)

sample.snpset <- new("oligoSnpSet", copyNumber=assayData(sample.snpset)[["copyNumber"]],
	   cnConfidence=assayData(sample.snpset)[["cnConfidence"]],
	   call=assayData(sample.snpset)[["calls"]],
	   callProbability=assayData(sample.snpset)[["callsConfidence"]],
	   annotation=annotation(sample.snpset),
	   phenoData=phenoData(sample.snpset),
	   featureData=featureData(sample.snpset),
	   experimentData=experimentData(sample.snpset))

##fD <- featureData(sample.snpset)
tmp <- new("AnnotatedDataFrame",
           data=data.frame(row.names=featureNames(sample.snpset)),
           varMetadata=data.frame(labelDescription=character()),
           dimLabels=c("featureNames", "featureColumns"))


tmp <- new("AnnotatedDataFrame")
tmp <- new("AnnotatedDataFrame",
           data=data.frame(),
           varMetadata=data.frame(),
           dimLabels=c("featureNames", "featureColumns"),
           
           
