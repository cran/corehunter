#' @import rJava
.onLoad <- function(libname, pkgname) {
  # init rJava
  .jpackage(pkgname, lib.loc = libname)
  # check Java version
  req.version <- 8
  version <- java.version()
  if(version < req.version){
    stop(sprintf(
      "Java version %d or later required. Found version %d.",
      req.version, version
    ))
  }
}
