.First.lib <- function(lib, pkgname) {
#  require(matchprobes)
  
  ## load the data
  thepath = system.file(package=pkgname)
  where   =  as.environment(match(paste("package:", pkgname, sep = ""),search()))

  data(list = pkgname, package = pkgname, envir=where)
}
