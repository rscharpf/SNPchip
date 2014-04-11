THISPKG <- "SNPchip"
.onAttach <- function(libname, pkgname) {
	packageStartupMessage("Welcome to SNPchip version ", packageDescription(THISPKG, fields="Version"))
}
.SNPchipPkgEnv <- new.env(parent=emptyenv())
