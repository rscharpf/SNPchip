setMethod("alleleA", "AnnotatedDataFrame", function(object) pData(object)$allele_a)
setMethod("alleleB", "AnnotatedDataFrame", function(object) pData(object)$allele_b)
##setMethod("chromosome", "AnnotatedDataFrame", function(object) pData(object)$chrom)
setMethod("dbSnpId", "AnnotatedDataFrame", function(object) pData(object)$dbsnp_rs_id)
##setMethod("position", "AnnotatedDataFrame", function(object) pData(object)$physical_pos)
##setMethod("enzyme", "AnnotatedDataFrame", function(object) pData(object)$enzyme)
setMethod("fragmentLength", "AnnotatedDataFrame", function(object) pData(object)$fragment_length)

