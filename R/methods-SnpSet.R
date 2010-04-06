setMethod("summary", signature(object = "SnpSet"),
          function(object, digits = 3, noCalls = FALSE, ...){
            het <- colMeans(ifelse(calls(object) == 2, 1, 0))
            hom <- colMeans(ifelse(calls(object) == 1 | calls(object) == 3, 1, 0))
            mat <- rbind(het, hom)
            mat
          })
