#### helper tools for export

#' @title Compute z-score corresponding to power
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
#' @examples # mean powers corresponding to the mean components of EM7z and 19.9
#' z_to_power(0:6, alpha = .05)
z_to_power  <- function(z, alpha = .05, a = stats::qnorm(alpha/2,lower.tail = FALSE)){
  if(!all(sapply(z, function(x)x >= 0)))stop("z must be >= 0")
  if(a  < 0)stop("a must be >= 0")
  if(is.null(a) & is.null(alpha))stop("Either 'alpha' or 'a' must be provided")
  if(is.null(alpha) & !is.null(a))alpha <- stats::pnorm(a, lower.tail = FALSE)*2
  if(alpha < 0 | alpha > 1)stop("alpha must be >= 0 & <= 1")
  1 - stats::pnorm(a, z, 1) + stats::pnorm(-a, z, 1)
}


#' @title Compute power corresponding to z-scores
#' @description A function for computing z-scores of two-sided tests
#' corresponding to power for a given significance level.
#' @param power A vector of z-scores
#' @param alpha Level of significance alpha
#' @param a Or, alternatively a z-score corresponding to \code{alpha}
#' @param nleqslv_control A named list of control parameters passed to the 
#' \link[nleqslv]{nleqslv} function used for solving the inverse of 
#' \link[=z_to_power]{z_to_power} function.  
#'
#' @export power_to_z
#'
#' @examples # z-scores corresponding to the (aproximate) power of components of EM7p
#' power_to_z(c(0.05, 0.20, 0.40, 0.60, 0.80, 0.974, 0.999), alpha = .05)
power_to_z  <- function(power, alpha = .05, a = stats::qnorm(alpha/2,lower.tail = FALSE),
                        nleqslv_control = list(xtol = 1e-15, maxit = 300, stepmax = .5)){
  if(!all(sapply(power, function(x)x >= 0 & x <= 1)))stop("power must be >= alpha & <= 1")
  if(a  < 0)stop("a must be >= 0")
  if(is.null(a) & is.null(alpha))stop("Either 'alpha' or 'a' must be provided")
  if(is.null(alpha) & !is.null(a))alpha <- stats::pnorm(a, lower.tail = FALSE)*2
  if(alpha < 0 | alpha > 1)stop("alpha must be >= 0 & <= 1")
  sapply(power, function(pow)nleqslv::nleqslv(.5, .solve_power_to_z, power = pow, a = a, control = nleqslv_control)$x)
}

.solve_power_to_z <- function(x, power, a){
  y = numeric(1)
  y = z_to_power(z = x, a = a) - power
  y
}

### internal tools for the results computation
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
  N_obs - .get_expected_N(EDR, N_sig)
}


# rounding for plot
.r2d <- function(x)format(round(x, 2), nsmall = 2)
.rXd <- function(x,X)format(round(x, X), nsmall = X)
