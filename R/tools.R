#### helper tools for export

#' @title Compute power corresponding to z-scores
#' @description A function for computing power of two-sided tests
#' corresponding to z-scores for a given significance level.
#' \code{alpha} (or corresponding cut-off z-score \code{a})
#'
#' @param z A vector of z-scores
#' @param alpha Level of significance alpha
#' @param a Or, alternatively a z-score corresponding to \code{alpha}
#'
#' @export z_to_power
#'
#' @examples # mean powers corresponding to the mean components of KD2
#' z_to_power(0:6, alpha = .05)
z_to_power  <- function(z, alpha = .05, a = stats::qnorm(alpha/2,lower.tail = FALSE)){
  if(!all(sapply(z, function(x)x >= 0)))stop("z must be >= 0")
  if(a  < 0)stop("a must be >= 0")
  if(is.null(a) & is.null(alpha))stop("Either 'alpha' or 'a' must be provided")
  if(is.null(alpha) & !is.null(a))alpha <- stats::pnorm(a, lower.tail = FALSE)*2
  if(alpha < 0 | alpha > 1)stop("alpha must be >= 0 & <= 1")
  1 - stats::pnorm(a, z, 1) + stats::pnorm(-a, z, 1)
}


#' @title Compute z-score corresponding to a power
#' @description A function for computing z-scores of two-sided tests
#' corresponding to power \code{power} for a given significance level 
#' alpha \code{alpha} (or corresponding cut-off z-statistic \code{a}).
#' 
#' @param power A vector of powers
#' @param alpha Level of significance alpha
#' @param a Or, alternatively a z-score corresponding to \code{alpha}
#' @param nleqslv_control A named list of control parameters passed to the 
#' \link[nleqslv]{nleqslv} function used for solving the inverse of 
#' \link[=z_to_power]{z_to_power} function.  
#'
#' @export power_to_z
#'
#' @examples # z-scores corresponding to the (aproximate) power of components of EM2
#' power_to_z(c(0.05, 0.20, 0.40, 0.60, 0.80, 0.974, 0.999), alpha = .05)
power_to_z  <- function(power, alpha = .05, a = stats::qnorm(alpha/2,lower.tail = FALSE),
                        nleqslv_control = list(xtol = 1e-15, maxit = 300, stepmax = .5)){
  if(a  < 0)stop("a must be >= 0")
  if(is.null(a) & is.null(alpha))stop("Either 'alpha' or 'a' must be provided")
  if(is.null(alpha) & !is.null(a))alpha <- stats::pnorm(a, lower.tail = FALSE)*2
  if(alpha < 0 | alpha > 1)stop("alpha must be >= 0 & <= 1")
  if(!all(sapply(power, function(x)x >= alpha & x <= 1)))stop("power must be >= alpha & <= 1")
  sapply(power, function(pow)nleqslv::nleqslv(.5, .solve_power_to_z, power = pow, a = a, control = nleqslv_control)$x)
}

.solve_power_to_z <- function(x, power, a){
  y = numeric(1)
  y = z_to_power(z = x, a = a) - power
  y
}

### internal tools for the results computation
.p_to_z     <- function(p){
  if(!all(sapply(p, function(x)x >= 0 & x <= 1)))stop("p-values must be >= 0 & <= 1")
  stats::qnorm(p/2, lower.tail = F)
}
# expected number of unpublished studies
.get_E_null <- function(N_fit, weights, Pow){
  sum((N_fit*weights)*(1/Pow-1))
}
# observed number of p-values > b
.get_O_high <- function(N_sig, N_fit){
  N_sig - N_fit
}
.get_EDR    <- function(N_sig, E_null){
  (N_sig)/(E_null + N_sig)
}
.get_ERR    <- function(N_fit, weights, Pow, O_high, N_sig){
  (sum(N_fit*weights*Pow) + O_high)/N_sig
}
.get_Z0     <- function(weights, N_fit, N_sig){
  weights[1]*(N_fit/N_sig)
}

.get_estimates <- function(z, a, N_fit, mu, weights){
  N_sig  <- length(z[z > a])

  Pow    <- z_to_power(z = mu, a = a)

  E_null <- .get_E_null(N_fit = N_fit, weights = weights, Pow = Pow)
  O_high <- .get_O_high(N_sig = N_sig, N_fit = N_fit)

  EDR <- .get_EDR(N_sig = N_sig, E_null = E_null)
  ERR <- .get_ERR(N_fit = N_fit, weights = weights, Pow = Pow, O_high = O_high, N_sig = N_sig)
  Z0  <- .get_Z0(weights = weights, N_fit = N_fit, N_sig = N_sig)

  estimates <- c(
    "ERR" = ERR,
    "EDR" = EDR,
    "Z0"  = Z0
  )
  return(estimates)
}

# additional functions for summary computation
.get_Soric_FDR     <- function(EDR, sig_level){
  ((1/EDR) - 1)*(sig_level/(1-sig_level))
}
.get_file_drawer_R <- function(EDR){
 (1-EDR)/EDR
}
.get_expected_N    <- function(EDR, N_sig){
  .get_file_drawer_R(EDR)*N_sig + N_sig
}
.get_missing_N     <- function(EDR, N_sig, N_obs){
  .get_expected_N(EDR, N_sig) - N_obs
}


# rounding for plot
.r2d  <- function(x)format(round(x, 2), nsmall = 2)
.rXd  <- function(x,X)format(round(x, X), nsmall = X)
.rXdn <- function(x,X)as.numeric(.rXd(x,X))

### the exported function

#' @title Reports whether x is a zcurve object
#'
#' @param x an object to test
#' @param ... additional arguments
#' @export is.zcurve
is.zcurve  <- function(x){
  inherits(x, "zcurve")
}

#' Prints estimates from z-curve object
#' @param x Estimate of a z-curve object
#' @param ... Additional arguments
#' @method print.estimates zcurve
#' @export print.summary.zcurve
#' @rawNamespace S3method(print, estimates.zcurve)
#' @seealso [zcurve()]
print.estimates.zcurve <- function(x){
  
  est_names  <- names(x[1:(length(x)-1)])
  est_values <- .rXdn(unlist(x[1:(length(x)-1)]), x$round.coef)
  names(est_values) <- est_names
  print(est_values)
  
}

#' @title z-curve estimates
#'
#' @description The following functions extract estimates 
#' from the z-curve object.
#' 
#' @param{object} the z-curve object
#' @param{round.coef} rounding for the printed values
#'
#' @export ERR
#' @export EDR
#' @export ODR
#' @export Soric
#' @export file_drawer_ration
#' @export expected_n
#' @export missing_n
#' @export significant_n
#' @export included_n
#' @name zcurve.estimates
#'
#' @details Technically, ODR, significant n, and included n 
#' are not z-curve estimates but they are grouped in this 
#' category for convenience.
#' @seealso [zcurve()]
NULL

#' @rdname zcurve.estimates
ERR   <- function(object, round.coef = 3){
  
  if(!is.zcurve(object))stop("The functions requires an 'zcurve' object.")

  sum <- summary(object)$coefficients

  val <- list()
  val[["Estimate"]] <- sum["ERR",1]
  
  if(!is.null(object[["boot"]])){
    val[["l.CI"]] <- sum["ERR", "l.CI"]
    val[["u.CI"]] <- sum["ERR", "u.CI"]
  }

  val[["round.coef"]] <- round.coef
  
  class(val) <- "estimates.zcurve"
  return(val)
}
#' @rdname zcurve.estimates
EDR   <- function(object, round.coef = 3){
  
  if(!is.zcurve(object))stop("The functions requires an 'zcurve' object.")
  
  sum <- summary(object)$coefficients
  
  val <- list()
  val[["Estimate"]] <- sum["EDR",1]
  
  if(!is.null(object[["boot"]])){
    val[["l.CI"]] <- sum["EDR", "l.CI"]
    val[["u.CI"]] <- sum["EDR", "u.CI"]
  }
  
  val[["round.coef"]] <- round.coef
  
  class(val) <- "estimates.zcurve"
  return(val)
}
#' @rdname zcurve.estimates
ODR   <- function(object, round.coef = 3){
  
  if(!is.zcurve(object))stop("The functions requires an 'zcurve' object.")
  
  sum <- summary(object)$model
  prt <- stats::prop.test(sum$N_sig, sum$N_all)
  
  val <- list()
  val[["Estimate"]] <- prt$estimate
  val[["l.CI"]] <- prt$conf.int[1]
  val[["u.CI"]] <- prt$conf.int[2]

  val[["round.coef"]] <- round.coef
  
  class(val) <- "estimates.zcurve"
  return(val)
}
#' @rdname zcurve.estimates
Soric <- function(object, round.coef = 3){
  
  if(!is.zcurve(object))stop("The functions requires an 'zcurve' object.")
  
  sum <- summary(object, all = TRUE)$coefficients
  
  val <- list()
  val[["Estimate"]] <- sum["Soric FDR",1]
  
  if(!is.null(object[["boot"]])){
    val[["l.CI"]] <- sum["Soric FDR", "l.CI"]
    val[["u.CI"]] <- sum["Soric FDR", "u.CI"]
  }
  
  val[["round.coef"]] <- round.coef
  
  class(val) <- "estimates.zcurve"
  return(val)
}
#' @rdname zcurve.estimates
file_drawer_ration <- function(object, round.coef = 3){
  
  if(!is.zcurve(object))stop("The functions requires an 'zcurve' object.")
  
  sum <- summary(object, all = TRUE)$coefficients
  
  val <- list()
  val[["Estimate"]] <- sum["File Drawer R",1]
  
  if(!is.null(object[["boot"]])){
    val[["l.CI"]] <- sum["File Drawer R", "l.CI"]
    val[["u.CI"]] <- sum["File Drawer R", "u.CI"]
  }
  
  val[["round.coef"]] <- round.coef
  
  class(val) <- "estimates.zcurve"
  return(val)
}
#' @rdname zcurve.estimates
expected_n         <- function(object, round.coef = 0){
  
  if(!is.zcurve(object))stop("The functions requires an 'zcurve' object.")
  
  sum <- summary(object, all = TRUE)$coefficients
  
  val <- list()
  val[["Estimate"]] <- sum["Expected N",1]
  
  if(!is.null(object[["boot"]])){
    val[["l.CI"]] <- sum["Expected N", "l.CI"]
    val[["u.CI"]] <- sum["Expected N", "u.CI"]
  }
  
  val[["round.coef"]] <- round.coef
  
  class(val) <- "estimates.zcurve"
  return(val)
}
#' @rdname zcurve.estimates
missing_n          <- function(object, round.coef = 0){
  
  if(!is.zcurve(object))stop("The functions requires an 'zcurve' object.")
  
  sum <- summary(object, all = TRUE)$coefficients
  
  val <- list()
  val[["Estimate"]] <- sum["Missing N",1]
  
  if(!is.null(object[["boot"]])){
    val[["l.CI"]] <- sum["Missing N", "l.CI"]
    val[["u.CI"]] <- sum["Missing N", "u.CI"]
  }
  
  val[["round.coef"]] <- round.coef
  
  class(val) <- "estimates.zcurve"
  return(val)
}
#' @rdname zcurve.estimates
significant_n      <- function(object){
  
  if(!is.zcurve(object))stop("The functions requires an 'zcurve' object.")
  
  sum <- summary(object)$model
  
  val <- list()
  val[["N"]] <- sum$N_sig
  
  val[["round.coef"]] <- 0
  
  class(val) <- "estimates.zcurve"
  return(val)
}
#' @rdname zcurve.estimates
included_n         <- function(object){
  
  if(!is.zcurve(object))stop("The functions requires an 'zcurve' object.")
  
  sum <- summary(object)$model
  
  val <- list()
  val[["N"]] <- sum$N_used

  val[["round.coef"]] <- 0
  
  class(val) <- "estimates.zcurve"
  return(val)
}
