NAMESPACE <- environment()
.module <- new("Module")

.onLoad <- function(libname, pkgname) {
  require(methods)
  unlockBinding(".module", NAMESPACE)
  assign(".module", Rcpp:::Module("rtm"), NAMESPACE)
  lockBinding(".module", NAMESPACE)
}

