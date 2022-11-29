#' @title Options for the zcurve package
#'
#' @description A placeholder object and functions for the zcurve package.
#' (adapted from the runjags R package).
#'
#' @param name the name of the option to get the current value of - for a list of
#' available options, see details below.
#' @param ... named option(s) to change - for a list of available options, see
#' details below.
#'
#' @return The current value of all available zcurve options (after applying any
#' changes specified) is returned invisibly as a named list.
#'
#' @export zcurve.options
#' @export zcurve.get_option
#' @name zcurve_options
#' @aliases zcurve_options zcurve.options zcurve.get_option
NULL


#' @rdname zcurve_options
zcurve.options    <- function(...){
  
  opts <- list(...)
  
  for(i in seq_along(opts)){
    
    if(!names(opts)[i] %in% names(zcurve.private))
      stop(paste("Unmatched or ambiguous option '", names(opts)[i], "'", sep=""))
    
    assign(names(opts)[i], opts[[i]] , envir = zcurve.private)
  }
  
  return(invisible(zcurve.private$options))
}

#' @rdname zcurve_options
zcurve.get_option <- function(name){
  
  if(length(name)!=1)
    stop("Only 1 option can be retrieved at a time")
  
  if(!name %in% names(zcurve.private))
    stop(paste("Unmatched or ambiguous option '", name, "'", sep=""))
  
  # Use eval as some defaults are put in using 'expression' to avoid evaluating at load time:
  return(eval(zcurve.private[[name]]))
}


zcurve.private <- new.env()
# Use 'expression' for functions to avoid having to evaluate before the package is fully loaded:
assign("defaultoptions",  list(envir    = zcurve.private))

assign("options",         zcurve.private$defaultoptions,   envir = zcurve.private)
assign("max_cores",       parallel::detectCores(logical = TRUE) - 1,  envir = zcurve.private)
