setMethod("unsplitSnpSet", c("list", "AnnotatedDataFrame"),
          function(from, annotatedDataFrame, ...){
            object <- switch(class(from[[1]]),
                             oligoSnpSet=new("oligoSnpSet", ...),
                             SnpCallSet=new("SnpCallSet", ...),
                             SnpCopyNumberSet=new("SnpCopyNumberSet",  ...),
                             stop("Not a defined class"))
            featureData(object) <- annotatedDataFrame[match(featureNames(object), featureNames(annotatedDataFrame)), ]
            phenoData(object) <- phenoData(from[[1]])
            annotation(object) <- annotation(from[[1]])
            experimentData(object) <- experimentData(from[[1]])
            stopifnot(validObject(object))
            object
          })


