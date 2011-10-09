require(Rcpp)

.onLoad <- function(libname, pkgname) {
  require(methods)
  loadRcppModules()
}

