#' @title Fit a z-curve to clustered data
#' 
#' @description \code{zcurve_clustered} is used to fit z-curve models to 
#' clustered data. The function requires a data object created with the 
#' [zcurve_data()] function as the input (where id denotes clusters). 
#' Two different methods that account for clustering ar implemented via 
#' the EM model: \code{"w"} for down weighting the likelihood of the test 
#' statistics proportionately to the number of repetitions in the clusters, 
#' and \code{"b"} for a nested bootstrap where only a single study from each 
#' bootstrap is selected for model fitting.
#' 
#' @param data an object created with [zcurve_data()] function.
#' @param method the method to be used for fitting. Possible options are
#' down weighting \code{"w"} and nested bootstrap \code{"b"}. 
#' Defaults to \code{"w"}.
#' @param bootstrap the number of bootstraps for estimating CI. To skip
#' bootstrap specify \code{FALSE}.
#' @param control additional options for the fitting algorithm more details in
#' \link[=control_EM]{control EM}.
#'
#'  
#' @references
#' \insertAllCited{}
#'
#' @return The fitted z-curve object
#'
#' @seealso [zcurve()], [summary.zcurve()], [plot.zcurve()], [control_EM], [control_density]
#' @export
zcurve_clustered <- function(data, method = "w", bootstrap = 1000, control = NULL){
  
  if(!method %in% c("w", "b"))
    stop("Wrong method, select a supported option.")
  if(method == "b" && is.logical(bootstrap) && !bootstrap)
    stop("The boostrap method requires bootstrap.")
  
  # set bootstrap
  if(!is.numeric(bootstrap)){
    bootstrap <- FALSE
  }else if(bootstrap <= 0){
    bootstrap <- FALSE
  }
  
  if(!inherits(data, "zcurve_data")){
    stop("The 'data' input must be created by the `zcurve_data()` function. See `?zcurve_data()` for more information.")
  }
  
  # create results object
  object            <- NULL
  object$call       <- match.call()
  object$method     <- switch(method, "w" = "EM (weighted)", "b" = "EM (bootstrapped)")
  object$input_type <- "zcurve-data"
  
  # create control
  control <- .zcurve_EM.control(control)
  
  ### prepare data
  if(nrow(data$precise) != 0){
    z    <- .p_to_z(data$precise$p)
    z_id <- data$precise$id
  }else{
    z    <- numeric()
    z_id <- numeric()
  }
  
  if(nrow(data$censored) != 0){
    
    lb   <- .p_to_z(data$censored$p.ub)
    ub   <- .p_to_z(data$censored$p.lb)
    b_id <- data$censored$id
    
    # remove non-significant censored p-values
    if(any(lb < control$a)){
      warning(paste0(sum(lb < control$a), " censored p-values removed due to the upper bound being larger that the fitting range."), immediate. = TRUE, call. = FALSE)
      b_id <- b_id[lb >= control$a]
      ub   <- ub[lb >= control$a]
      lb   <- lb[lb >= control$a]
    }
    
    # move too significant censored p-values among precise p-values          
    if(length(lb) > 0 && any(lb >= control$b)){
      object$data <- c(object$data, lb[lb >= control$b])
      b_id <- b_id[lb < control$b]
      ub   <- ub[lb < control$b]
      lb   <- lb[lb < control$b]
    }
    
    if(length(lb) > 0){
      # restrict the upper censoring to the fitting range 
      ub <- ifelse(ub > control$b, control$b, ub)
  
      # update control
      control$type <- 3
    }else{
      lb   <- numeric()    
      ub   <- numeric()
      b_id <- numeric()
    }
  }else{
    lb   <- numeric()    
    ub   <- numeric()
    b_id <- numeric()
  }
  
  object$data           <- z
  object$data_id        <- z_id
  object$data_censoring <- data.frame(lb = lb, ub = ub, id = b_id)
  
  object$control        <- control
  
  # only run the algorithm with some significant results
  if(sum(z > control$a & z < control$b) + length(lb) < 10)
    stop("There must be at least 10 z-scores in the fitting range but a much larger number is recommended.")
  
  
  # use appropriate algorithm
  if(method == "b"){
    fit_b <- .zcurve_EM_b(z = z, z_id = z_id, lb = lb, ub = ub, b_id = b_id, control = control, bootstrap = bootstrap)
    fit   <- fit_b$fit
  }else if(method == "w"){
    fit   <- .zcurve_EM_w(z = z, z_id = z_id, lb = lb, ub = ub, b_id = b_id, control = control)
  }
  object$fit <- fit
  
  # check convergence
  object$converged <- ifelse(fit$iter < control$max_iter, TRUE, FALSE)
  
  
  # do bootstrap
  if(bootstrap != FALSE){
    if(method == "b"){
      fit_boot <- fit_b$fit_boot
    }else if(method == "w"){
      fit_boot <- .zcurve_EM_w_boot(z = z, z_id = z_id, lb = lb, ub = ub, b_id = b_id, control = control, fit = fit, bootstrap = bootstrap) 
    }
    object$boot <- fit_boot
  }
  
  
  # estimates
  object$coefficients <- .get_estimates(mu = fit$mu, weights = fit$weights, prop_high = fit$prop_high, sig_level = control$sig_level, a = control$a)
  
  
  # boot estimates
  if(bootstrap != FALSE){
    object$coefficients_boot <- data.frame(t(sapply(1:bootstrap, function(i){
      .get_estimates(mu = fit_boot$mu[i,], weights = fit_boot$weights[i,], prop_high = fit_boot$prop_high[i], sig_level = control$sig_level, a = control$a)
    })))
  }
  
  
  class(object) <- c("zcurve", "zcurve.clustered")
  return(object)
}

.zcurve_EM_b          <- function(z, z_id, lb, ub, b_id, control, bootstrap){
  
  # get starting value z-curves
  fit_start <- .zcurve_EMc_start_fast_RCpp(x           = z,
                                           lb          = lb,
                                           ub          = ub,
                                           K           = control$K,
                                           mu          = control$mu,
                                           sigma       = control$sigma,
                                           mu_alpha    = control$mu_alpha,
                                           mu_max      = control$mu_max,
                                           theta_alpha = control$theta_alpha,
                                           a           = control$a,
                                           b           = control$b,
                                           sig_level   = control$sig_level,
                                           fit_reps    = control$fit_reps,
                                           max_iter    = control$max_iter_start,
                                           criterion   = control$criterion_start)
  
  # fit final z-curve
  data_index <- rbind(
    data.frame(
      z    = z,
      lb   = rep(NA, length(z)),
      ub   = rep(NA, length(z)),
      id   = z_id,
      type = rep(1, length(z))
    ),
    data.frame(
      z    = rep(NA, length(lb)),
      lb   = lb,
      ub   = ub,
      id   = b_id,
      type = rep(2, length(lb))
    )
  )
  
  
  fit <- list()
  for(i in 1:bootstrap){
    
    boot_data <- .boot_id(data_index)
    
    fit[[i]]  <- .zcurve_EMc_fit_fast_RCpp(x          = boot_data$z[boot_data$type == 1],
                                           lb         = boot_data$lb[boot_data$type == 2],
                                           ub         = boot_data$ub[boot_data$type == 2],
                                           mu         = fit_start$mu[which.max(fit_start$Q),],
                                           sigma      = control$sigma,
                                           theta      = fit_start$weights[which.max(fit_start$Q),],
                                           a          = control$a,
                                           b          = control$b,
                                           sig_level  = control$sig_level,
                                           max_iter   = control$max_iter,
                                           criterion  = control$criterion)
    
  }
  
  fit_boot = list(
    "mu"         = do.call(rbind, lapply(fit, function(f) f[["mu"]])),
    "weights"    = do.call(rbind, lapply(fit, function(f) f[["weights"]])),
    "prop_high"  = do.call(c, lapply(fit, function(f) f[["prop_high"]])),
    "Q"          = do.call(c, lapply(fit, function(f) f[["Q"]])),
    "iter"       = do.call(c, lapply(fit, function(f) f[["iter"]]))
  )
  
  return(
    list(
      fit =      list(
        "mu"         = apply(fit_boot[["mu"]], 2, mean),
        "weights"    = apply(fit_boot[["weights"]], 2, mean),
        "prop_high"  = mean(fit_boot[["prop_high"]]),
        "Q"          = mean(fit_boot[["Q"]]),
        "iter"       = mean(fit_boot[["iter"]]),
        "iter_start" = fit_start$iter[which.max(fit_start$Q)]
      ),
      fit_boot = fit_boot
    )
  )
}
.zcurve_EM_w          <- function(z, z_id, lb, ub, b_id, control){
  
  # get starting value z-curves
  fit_start <- .zcurve_EMc_start_fast_RCpp(x           = z,
                                           lb          = lb,
                                           ub          = ub,
                                           K           = control$K,
                                           mu          = control$mu,
                                           sigma       = control$sigma,
                                           mu_alpha    = control$mu_alpha,
                                           mu_max      = control$mu_max,
                                           theta_alpha = control$theta_alpha,
                                           a           = control$a,
                                           b           = control$b,
                                           sig_level   = control$sig_level,
                                           fit_reps    = control$fit_reps,
                                           max_iter    = control$max_iter_start,
                                           criterion   = control$criterion_start)
  
  # compute weights
  id_freq    <- table(c(z_id, b_id))
  id_weights <- data.frame(
    id = names(id_freq),
    w  = 1/as.vector(id_freq)
  ) 
  
  # fit final z-curve
  fit  <- .zcurve_EMc_fit_fast_w_RCpp(x          = z,
                                      x_w        = id_weights$w[match(z_id, id_weights$id)],
                                      lb         = lb,
                                      ub         = ub,
                                      b_w        = id_weights$w[match(b_id, id_weights$id)],
                                      mu         = fit_start$mu[which.max(fit_start$Q),],
                                      sigma      = control$sigma,
                                      theta      = fit_start$weights[which.max(fit_start$Q),],
                                      a          = control$a,
                                      b          = control$b,
                                      sig_level  = control$sig_level,
                                      max_iter   = control$max_iter,
                                      criterion  = control$criterion)
  
  return(
    list(
      "mu"         = fit$mu,
      "weights"    = fit$weights,
      "prop_high"  = fit$prop_high,
      "Q"          = fit$Q,
      "iter"       = fit$iter,
      "iter_start" = fit_start$iter[which.max(fit_start$Q)]
    )
  )
}
.zcurve_EM_w_boot     <- function(z, z_id, lb, ub, b_id, control, fit, bootstrap){
  
  # compute weights
  id_freq    <- table(c(z_id, b_id))
  id_weights <- data.frame(
    id = names(id_freq),
    w  = 1/as.vector(id_freq)
  ) 
  
  x_w <- id_weights$w[match(z_id, id_weights$id)]
  b_w <- id_weights$w[match(b_id, id_weights$id)]
  
  indx <- c(
    if(length(z)  > 0) 1:length(z),
    if(length(lb) > 0) (-length(lb)):-1
  )
  
  fit_boot <- .zcurve_EMc_boot_fast_w_RCpp(x         = z,
                                           x_w       = x_w,
                                           lb        = lb,
                                           ub        = ub,
                                           b_w       = b_w,
                                           indx      = indx,
                                           mu        = fit$mu,
                                           sigma     = control$sigma,
                                           theta     = fit$weights,
                                           a         = control$a,
                                           b         = control$b,
                                           sig_level = control$sig_level,
                                           bootstrap = bootstrap,
                                           criterion = control$criterion_boot,
                                           max_iter  = control$max_iter_boot)
  return(
    list(
      "mu"        = fit_boot$mu,
      "weights"   = fit_boot$weights,
      "Q"         = fit_boot$Q,
      "prop_high" = fit_boot$prop_high,
      "iter"      = fit_boot$iter
    )
  )
  
}

.boot_id <- function(data){
  
  unique_id <- data$id
  
  if(length(unique_id) == nrow(data)){
    return(data)
  }
  
  boot_out <- list()
  for(id in unique_id){
    
    temp_data <- data[data$id == id,]
    
    if(nrow(data) == 1){
      boot_out[[id]] <- temp_data
    }else{
      boot_out[[id]] <- temp_data[sample(1, nrow(temp_data)),]
    }
  }
  boot_out <- do.call(rbind, boot_out)
  
  return(boot_out)
}
