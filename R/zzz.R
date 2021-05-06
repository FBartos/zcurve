.onAttach <- function(libname, pkgname){
  
  packageStartupMessage(
    "Please, note the following changes in version 1.0.9 (see NEWS for more details):\n- The ERR estimate now takes the directionality of the expected replications into account, which might lead to slight changes in the estimates."
  )
}