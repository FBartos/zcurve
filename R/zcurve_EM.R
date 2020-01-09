.zcurve_EM          <- function(z_sig, control){

  # get starting value z-curves
  if(control$type == 1){
    fit_start <- .zcurve_EM_start_fast_RCpp(x           = z_sig,
                                            K           = control$K,
                                            mu          = control$mu,
                                            sigma       = control$sigma,
                                            mu_alpha    = control$mu_alpha,
                                            mu_max      = control$mu_max,
                                            theta_alpha = control$theta_alpha,
                                            a           = control$a,
                                            b           = control$b,
                                            fit_reps    = control$fit_reps,
                                            max_iter    = control$max_iter_start,
                                            criterion   = control$criterion_start)
  }else if(control$type == 2){
    fit_start <- .zcurve_EM_start_RCpp(x           = z_sig,
                                       type        = control$type,
                                       K           = control$K,
                                       mu          = control$mu,
                                       sigma       = control$sigma,
                                       mu_alpha    = control$mu_alpha,
                                       mu_max      = control$mu_max,
                                       theta_alpha = control$theta_alpha,
                                       a           = control$a,
                                       b           = control$b,
                                       fit_reps    = control$fit_reps,
                                       max_iter    = control$max_iter_start,
                                       criterion   = control$criterion_start)
  }

  # fit final z-curve
  if(control$type == 1){
    fit  <- .zcurve_EM_fit_fast_RCpp(x          = z_sig,
                                     mu         = fit_start$mu[which.max(fit_start$Q),],
                                     sigma      = control$sigma,
                                     theta      = fit_start$weights[which.max(fit_start$Q),],
                                     a          = control$a,
                                     b          = control$b,
                                     max_iter   = control$max_iter,
                                     criterion  = control$criterion)
  }else if(control$type == 2){
    fit  <- .zcurve_EM_fit_RCpp(x          = z_sig,
                                type       = control$type,
                                mu         = fit_start$mu[which.max(fit_start$Q),],
                                sigma      = control$sigma,
                                theta      = fit_start$weights[which.max(fit_start$Q),],
                                a          = control$a,
                                b          = control$b,
                                max_iter   = control$max_iter,
                                criterion  = control$criterion)
  }

  return(
    list(
      "mu"         = fit$mu,
      "weights"    = fit$weights,
      "N_fit"      = fit$N_fit,
      "Q"          = fit$Q,
      "iter"       = fit$iter,
      "iter_start" = fit_start$iter[which.max(fit_start$Q)]
    )
  )

}
.zcurve_EM_boot     <- function(z_sig, control, fit, bootstrap){

  if(control$type == 1){
    fit_boot <- .zcurve_EM_boot_fast_RCpp(x    = z_sig,
                                     mu        = fit$mu,
                                     sigma     = control$sigma,
                                     theta     = fit$weights,
                                     a         = control$a,
                                     b         = control$b,
                                     bootstrap = bootstrap,
                                     max_iter  = control$max_iter_boot,
                                     criterion = control$criterion_boot
                                     )
  }else if(control$type == 2){
    fit_boot <- .zcurve_EM_boot_RCpp(x         = z_sig,
                                     type      = control$type,
                                     mu        = fit$mu,
                                     sigma     = control$sigma,
                                     theta     = fit$weights,
                                     a         = control$a,
                                     b         = control$b,
                                     bootstrap = bootstrap,
                                     criterion = control$criterion_boot,
                                     max_iter  = control$max_iter_boot)
  }

  return(
    list(
      "mu"      = fit_boot$mu,
      "weights" = fit_boot$weights,
      "Q"       = fit_boot$Q,
      "N_fit"   = fit_boot$N_fit,
      "iter"    = fit_boot$iter
    )
  )

}

#' @name control_EM
#' @title Control settings for the zcurve EM algorithm
#' @description All these settings are passed to the Expectation Maximization
#' fitting algorithm. All unspecified settings are set to the default value.
#' Setting \code{model = "EM"} sets all settings to the default
#' value irrespective of any other setting and fits z-curve as described in
#' \insertCite{zcurve2;textual}{zcurve}
#' 
#' @param model A type of model to be fitted, defaults to \code{"EM"}
#' for a z-curve with 7 z-scores centered components.
#' @param sig_level An alpha level of the test statistics, defaults to
#' \code{.05}
#' @param a A beginning of fitting interval, defaults to
#' \code{qnorm(sig_level/2,lower.tail = F)}
#' @param b An end of fitting interval, defaults to \code{5}
#' @param mu Means of the components, defaults to
#' \code{0:6}
#' @param sigma A standard deviation of the components, defaults to
#' \code{rep(1, length(mu))}
#' @param theta_alpha A vector of alpha parameters of a Dirichlet distribution
#' for generating random starting values for the weights, defaults to
#' \code{rep(.5, length(mu))}
#' @param theta_max Upper limits for weights, defaults to
#' \code{rep(1,length(mu))}
#' @param criterion A criterion to terminate the EM algorithm,
#' defaults to \code{1e-6}
#' @param criterion_start A criterion to terminate the starting phase 
#' of the EM algorithm, defaults to \code{1e-3}
#' @param criterion_boot A criterion to terminate the bootstrapping phase 
#' of the EM algorithm, defaults to \code{1e-5}
#' @param max_iter A maximum number of iterations of the EM algorithm
#' (not including the starting iterations) defaults to \code{10000}
#' @param max_iter_start A maximum number of iterations for the 
#' starting phase of EM algorithm, defaults to \code{100}
#' @param max_iter_boot A maximum number of iterations for the 
#' booting phase of EM algorithm, defaults to \code{100}
#' @param fit_reps A number of starting fits to get the initial
#' position for the EM algorithm, defaults to \code{100}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @examples # to increase the number of starting fits
#' # and change the means of the mixture components
#' 
#' ctrl <- list(
#'    fit_reps  = 50,
#'    mu = c(0, 1.5, 3, 4.5, 6)
#' )
#' \dontrun{zcurve(z, algorithm = "EM", control = ctrl)}
#' 
#' @seealso [zcurve()], [control_density]
NULL

.zcurve_EM.control  <- function(control){
  #if(is.null(control)){
  #  control$sig_level       <- .05
  #  control$a               <- stats::qnorm(control$sig_level/2,lower.tail = F)
  #  control$b               <- 5
  #  control$type            <- 1 # legacy from z-curve with estimated means
  #  control$mu              <- c(0, 1.11, 1.71, 2.21, 2.80, 3.9, 5)
  #  control$sigma           <- rep(1, length(control$mu))
  #  control$K               <- length(control$mu) # legacy from z-curve with estimated means
  #  control$theta_alpha     <- rep(.5, length(control$mu))
  #  control$mu_alpha        <- 2:(control$K+1) # legacy from z-curve with estimated means
  #  control$mu_max          <- control$b + 2   # legacy from z-curve with estimated means
  #  control$criterion       <- 1e-5
  #  control$max_iter        <- 1000
  #  control$criterion_start <- 1e-3
  #  control$max_iter_start  <- 100
  #  control$fit_reps        <- 20
  #  control$model           <- "EM7p"
  #  return(control)
  #}
  if(is.null(control)){
    control$sig_level       <- .05
    control$a               <- stats::qnorm(control$sig_level/2,lower.tail = F)
    control$b               <- 6
    control$type            <- 1
    control$mu              <- 0:6
    control$sigma           <- rep(1, length(control$mu))
    control$K               <- length(control$mu)
    control$theta_alpha     <- rep(.5, length(control$mu))
    control$mu_alpha        <- 2:(control$K+1)
    control$mu_max          <- control$b + 2
    control$criterion       <- 1e-6
    control$max_iter        <- 10000
    control$criterion_boot  <- 1e-5
    control$max_iter_boot   <- 1000
    control$criterion_start <- 1e-3
    control$max_iter_start  <- 100
    control$fit_reps        <- 100
    control$model           <- "EM"
    return(control)
  }
  if(!is.null(control$model)){
    if(control$model == "EM"){
      control$sig_level       <- .05
      control$a               <- stats::qnorm(control$sig_level/2,lower.tail = F)
      control$b               <- 6
      control$type            <- 1
      control$mu              <- 0:6
      control$sigma           <- rep(1, length(control$mu))
      control$K               <- length(control$mu)
      control$theta_alpha     <- rep(.5, length(control$mu))
      control$mu_alpha        <- 2:(control$K+1)
      control$mu_max          <- control$b + 2
      control$criterion       <- 1e-6
      control$max_iter        <- 10000
      control$criterion_boot  <- 1e-5
      control$max_iter_boot   <- 1000
      control$criterion_start <- 1e-3
      control$max_iter_start  <- 100
      control$fit_reps        <- 100
      control$model           <- "EM"
      return(control)
    }
  }
  if(is.null(control$sig_level)){
    control$sig_level       <- .05
  }
  if(is.null(control$a)){
    control$a               <- stats::qnorm(control$sig_level/2,lower.tail = F)
  }
  if(is.null(control$b)){
    control$b               <- 6
  }
  if(is.null(control$type)){
    control$type            <- 1
  }
  if(is.null(control$mu)){
    control$mu              <- 0:6
  }
  if(is.null(control$sigma)){
    control$sigma           <- rep(1, length(control$mu))
  }
  if(is.null(control$K)){
    control$K               <- ifelse(control$type == 1, length(control$mu), 4)
  }
  if(is.null(control$theta_alpha)){
    control$theta_alpha     <- rep(.5, length(control$mu))
  }
  if(is.null(control$mu_alpha)){
    control$mu_alpha        <- 2:(control$K+1)
  }
  if(is.null(control$mu_max)){
    control$mu_max          <- control$b + 2
  }
  if(is.null(control$criterion)){
    control$criterion       <- 1e-6
  }
  if(is.null(control$max_iter)){
    control$max_iter        <- 10000
  }
  if(is.null(control$criterion_boot)){
    control$criterion_boot  <- 1e-5
  }
  if(is.null(control$max_iter_boot)){
    control$max_iter_boot   <- 1000
  }
  if(is.null(control$criterion_start)){
    control$criterion_start <- 1e-3
  }
  if(is.null(control$max_iter_start)){
    control$max_iter_start  <- 100
  }
  if(is.null(control$fit_reps)){
    control$fit_reps        <- 100
  }
  if(is.null(control$model)){
    control$model           <- NULL
  }
  return(control)
}
