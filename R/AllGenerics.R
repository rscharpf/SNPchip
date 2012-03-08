setGeneric("xyplot2", function(x, data, range, frame=50e3L, ...) standardGeneric("xyplot2"))
setGeneric("xyplot", useAsDefault=function(x, data, ...) lattice::xyplot(x, data, ...))
##setGeneric("cloud", useAsDefault=function(x, data, ...) lattice::cloud(x, data, ...))
