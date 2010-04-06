THISPKG <- "SNPchip"
.onLoad <- function(libname, pkgname) {
  require("methods")
}

.onAttach <- function(libname, pkgname) {
  message("Welcome to SNPchip version ", packageDescription(THISPKG, field="Version"))
}

##.onUnload <- function( libpath ){
##  library.dynam.unload(THISPKG, libpath)
##}

.SNPchipPkgEnv <- new.env(parent=emptyenv())
