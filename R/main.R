#' @title Fit a z-curve
#' 
#' @description \code{zcurve} is used to fit z-curve models. The function
#' takes input of z-statistics or two-sided p-values and returns object of
#' class \code{"zcurve"} that can be further interrogated by summary and plot
#' function. It default to EM model, but different version of z-curves can
#' be specified using the \code{method} and \code{control} arguments. See
#' 'Examples' and 'Details' for more information.
#' 
#' @param z a vector of z-scores.
#' @param p a vector of two-sided p-values, internally transformed to 
#' z-scores.
#' @param method the method to be used for fitting. Possible options are
#' Expectation Maximization \code{"EM"} and density \code{"density"},
#' defaults to \code{"EM"}.
#' @param bootstrap the number of bootstraps for estimating CI. To skip
#' bootstrap specify \code{FALSE}.
#' @param control additional options for the fitting algorithm more details in
#' \link[=control_EM]{control EM} or \link[=control_density]{control density}.
#'
#' @details The function returns the EM method by default and changing 
#' \code{method = "density"} gives the KD2 version of z-curve as outlined in
#' \insertCite{zcurve2;textual}{zcurve}. For the original z-curve 
#' \insertCite{zcurve1}{zcurve}, referred to as KD1, specify 
#'  \code{'control = "density", control = list(model = "KD1")'}.
#'  
#' @references
#' \insertAllCited{}
#'
#' @return The fitted z-curve object
#' @export zcurve
#'
#' @examples
#' # load data from OSC 2015 reproducibility project
#' OSC.z
#'
#' # fit an EM z-curve (with disabled bootstrap due to examples times limits)
#' m.EM <- zcurve(OSC.z, method = "EM", bootstrap = FALSE)
#' # a version with 1000 boostraped samples would looked like:
#' \donttest{m.EM <- zcurve(OSC.z, method = "EM", bootstrap = 1000)}
#' 
#' # or KD2 z-curve (use larger bootstrap for real inference)
#' m.D <- zcurve(OSC.z, method = "density", bootstrap = FALSE)
#'
#' # inspect the results
#' summary(m.EM)
#' summary(m.D)
#' # see '?summary.zcurve' for more output options
#' 
#' # plot the results
#' plot(m.EM)
#' plot(m.D)
#' # see '?plot.zcurve' for more plotting options
#'
#' # to specify more options, set the control arguments
#' # ei. increase the maximum number of iterations and change alpha level
#' ctr1 <- list(
#'   "max_iter" = 9999,
#'   "alpha"    = .10
#'   )
#' \donttest{m1.EM <- zcurve(OSC.z, method = "EM", bootstrap = FALSE, control = ctr1)}
#' # see '?control_EM' and '?control_density' for more information about different
#' # z-curves specifications
#' @seealso [summary.zcurve()], [plot.zcurve()], [control_EM], [control_density]
zcurve       <- function(z, p, method = "EM", bootstrap = 1000, control = NULL){
  
  # check input
  input_type <- NULL
  if(missing(z) & missing(p))stop("No data input")
  if(!missing(z)){
    if(!is.numeric(z))stop("Wrong z-scores input: Data are not nummeric.")
    if(!is.vector(z))stop("Wrong z-scores input: Data are not a vector")
    if(all(z <= 1 & z >= 0))stop("It looks like you are entering p-values rather than z-scores. To use p-values, explicitly name your argument 'zcurve(p = [vector of p-values])'")
    input_type <- c(input_type, "z")
  }else{
    z <- NULL
  }
  if(!missing(p)){
    if(!is.numeric(p))stop("Wrong p-values input: Data are not nummeric.")
    if(!is.vector(p))stop("Wrong p-values input: Data are not a vector") 
    input_type <- c(input_type, "p")
  }else{
    p <- NULL
  }
  
  if(!method %in% c("EM", "density"))stop("Wrong method, select a supported option")
  if(!is.numeric(bootstrap))bootstrap <- FALSE
  if(bootstrap <= 0)        bootstrap <- FALSE
  
  
  # create results object
  object            <- NULL
  object$call       <- match.call()
  object$method     <- method
  object$input_type <- input_type
  
  
  # update control
  if(method == "EM"){
    control <- .zcurve_EM.control(control)
  }else if(method == "density"){
    control <- .zcurve_density.control(control)
  }
  object$control <- control
  
  
  # prepare data
  if(!is.null(p)){
    z_from_p <- .p_to_z(p)
    z        <- c(z, z_from_p)
  }
  z           <- abs(z)
  object$data <- z
  
  # only run the algorithm with some significant results
  if(sum(z > control$a & z < control$b) < 10)stop("There must be at least 10 z-scores in the fitting range but a much larger number is recommended.")
  
  # use apropriate algorithm
  if(method == "EM"){
    fit <- .zcurve_EM(z = z, control = control)
  }else if(method == "density"){
    fit <- .zcurve_density(z = z, control = control)
  }
  object$fit <- fit
  
  
  # check convergence
  if(method == "EM"){
    object$converged <- ifelse(fit$iter < control$max_iter, TRUE, FALSE)
  }else if(method == "density"){
    object$converged <- fit$converged
    if(fit$message == "singular convergence (7)")object$converged <- TRUE
    if(fit$message == "both X-convergence and relative convergence (5)")object$converged <- TRUE
  }
  if(object$converged == FALSE)warning("Model did not converge.")
  
  
  # do bootstrap
  if(bootstrap != FALSE){
    # use apropriate algorithm
    if(method == "EM"){
      fit_boot <- .zcurve_EM_boot(z = z, control = control, fit = fit, bootstrap = bootstrap)
    }else if(method == "density"){
      fit_boot <- .zcurve_density_boot(z = z, control = control, bootstrap = bootstrap)
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
  
  
  class(object) <- "zcurve"
  return(object)
}

### methods
#' Prints a fitted z-curve object
#' @param x Fitted z-curve object
#' @param ... Additional arguments
#' @export  print.zcurve
#' @rawNamespace S3method(print, zcurve)
#' @seealso [zcurve()]
print.zcurve         <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nEstimates:\n")
  print(x$coefficients[1:2])
}

#' Summarize fitted z-curve object
#'
#' @param object A fitted z-curve object.
#' @param type Whether the results \code{"results"} or the 
#' mixture mode parameters \code{"parameters"} should be 
#' returned. Defaults to  \code{"results"}.
#' @param all Whether additional results, such as file drawer 
#' ration, expected and missing number of studies, and Soric FDR 
#' be returned. Defaults to \code{FALSE}
#' @param ERR.adj Confidence intervals adjustment for ERR. Defaults 
#' to \code{.03} as proposed by Bartos & Schimmack (in preparation).
#' @param EDR.adj Confidence intervals adjustment for EDR. Defaults 
#' to \code{.05} as proposed by Bartos & Schimmack (in preparation).
#' @param round.coef To how many decimals should the coefficient 
#' be rounded. Defaults to \code{3}.
#' @param ... Additional arguments
#'
#' @return Summary of a z-curve object
#' 
#' @method summary zcurve
#' @export summary.zcurve
#' @rawNamespace S3method(summary, zcurve)
#' @seealso [zcurve()]
summary.zcurve       <- function(object, type = "results", all = FALSE, ERR.adj = .03, EDR.adj = .05, round.coef = 3, ...){

  if(object$method == "EM"){
    if(!is.null(object$boot)){
      fit_index <- c(object$fit$Q, unname(stats::quantile(object$boot$Q, c(.025, .975))))
    }else{
      fit_index <- c(object$fit$Q)
    }
    method_text <- object$method
    iter_text   <- paste(c(object$fit$iter_start," + ",object$fit$iter), collapse = "")
    fit_stat    <- "Q"
  }else if(object$method == "density"){
    if(!is.null(object$boot)){
      fit_index <- c(object$fit$objective, unname(stats::quantile(object$boot$objective, c(.025, .975))))
    }else{
      fit_index <- c(object$fit$objective)
    }
    if(object$control$version == 1){
      method_text <- paste(c(object$method, " (version 1)"), collapse = "")
      fit_stat    <- "MAE (*1e3)"
      fit_index   <- fit_index*1e3
    }else{
      method_text <- object$method
      fit_stat    <- "RMSE"
    }
    iter_text <- object$fit$iter
  }
  
  temp_N_sig     <- sum(object$data > stats::qnorm(object$control$sig_level/2, lower.tail = FALSE))
  temp_N_obs     <- length(object$data)
  temp_N_used    <- sum(object$data > object$control$a & object$data < object$control$b)
    
  model <- list(
    "method"    = method_text,
    "model"     = ifelse(is.null(object$control$model), "custom", object$control$model),
    "fit_index" = fit_index,
    "fit_stat"  = fit_stat,
    "iter"      = iter_text,
    "input_type"= object$input_type,
    "N_all"     = temp_N_obs,
    "N_sig"     = temp_N_sig,
    "N_used"    = temp_N_used
  )
  
  if(type == "results" | substr(type,1,3) == "res"){
    
    if(!is.null(object$boot)){
      l.CI <- c(stats::quantile(object$coefficients_boot$ERR, .025),
                stats::quantile(object$coefficients_boot$EDR, .025))
      u.CI <- c(stats::quantile(object$coefficients_boot$ERR, .975),
                stats::quantile(object$coefficients_boot$EDR, .975))
    }else{
      l.CI <- NULL
      u.CI <- NULL
    }
    
    TAB <- cbind(Estimate = stats::coef(object)[1:2],
                 l.CI     = l.CI,
                 u.CI     = u.CI)
    
    # adjust CIs
    if(!is.null(object$boot)){
      TAB["ERR", "l.CI"] <- ifelse(TAB["ERR", "l.CI"] - ERR.adj < object$control$sig_level, object$control$sig_level, TAB["ERR", "l.CI"] - ERR.adj)
      TAB["ERR", "u.CI"] <- ifelse(TAB["ERR", "u.CI"] + ERR.adj > 1, 1, TAB["ERR", "u.CI"] + ERR.adj)
      TAB["EDR", "l.CI"] <- ifelse(TAB["EDR", "l.CI"] - EDR.adj < object$control$sig_level, object$control$sig_level, TAB["EDR", "l.CI"] - EDR.adj)
      TAB["EDR", "u.CI"] <- ifelse(TAB["EDR", "u.CI"] + EDR.adj > 1, 1, TAB["EDR", "u.CI"] + EDR.adj)
      TAB["EDR", "u.CI"] <- ifelse(TAB["EDR", "u.CI"] > TAB["ERR", "u.CI"], TAB["ERR", "u.CI"], TAB["EDR", "u.CI"])
    }
    
    # additional stats
    if(all){
      
      temp_sig_level <- object$control$sig_level
      
      if(!is.null(object$boot)){
        TAB <- rbind(
          TAB,
          "Soric FDR"     = c(
            "Estimate"     = .get_Soric_FDR(TAB["EDR","Estimate"], temp_sig_level),
            "l.CI"         = .get_Soric_FDR(TAB["EDR","u.CI"],     temp_sig_level),
            "u.CI"         = .get_Soric_FDR(TAB["EDR","l.CI"],     temp_sig_level)),
          "File Drawer R" = c(
            "Estimate"     = .get_file_drawer_R(TAB["EDR","Estimate"]),
            "l.CI"         = .get_file_drawer_R(TAB["EDR","u.CI"]),
            "u.CI"         = .get_file_drawer_R(TAB["EDR","l.CI"])),
          "Expected N"    = c(
            "Estimate"     = .get_expected_N(TAB["EDR","Estimate"], temp_N_sig),
            "l.CI"         = .get_expected_N(TAB["EDR","u.CI"],     temp_N_sig),
            "u.CI"         = .get_expected_N(TAB["EDR","l.CI"],     temp_N_sig)),
          "Missing N"     = c(
            "Estimate"     = .get_missing_N(TAB["EDR","Estimate"], temp_N_sig, temp_N_obs),
            "l.CI"         = .get_missing_N(TAB["EDR","u.CI"],     temp_N_sig, temp_N_obs),
            "u.CI"         = .get_missing_N(TAB["EDR","l.CI"],     temp_N_sig, temp_N_obs))
        )
      }else{
        TAB <- rbind(
          TAB,
          "Soric FDR"     = c("Estimate" = .get_Soric_FDR(TAB["EDR","Estimate"], temp_sig_level)),
          "File Drawer R" = c("Estimate" = .get_file_drawer_R(TAB["EDR","Estimate"])),
          "Expected N"    = c("Estimate" = .get_expected_N(TAB["EDR","Estimate"], temp_N_sig)),
          "Missing N"     = c("Estimate" = .get_missing_N(TAB["EDR","Estimate"], temp_N_sig, temp_N_obs)))
      }
    }
    
    if(object$method == "density"){
      if(object$control$version == 1){
        TAB <- as.data.frame(TAB[1,,drop=FALSE])
      }
    }
    
    
  }else if(type == "parameters" | substr(type,1,3) == "par"){
    
    if(!is.null(object$boot)){
      
      if(object$method == "density"){
        if(object$control$version == 1){
          M_l.CI <- apply(object$boot$mu,2,stats::quantile, prob = .025)
          M_u.CI <- apply(object$boot$mu,2,stats::quantile, prob = .975)
        }else{
          M_l.CI <- NULL
          M_u.CI <- NULL
        }
      }else{
        M_l.CI <- NULL
        M_u.CI <- NULL
      }

      W_l.CI <- apply(object$boot$weights,2,stats::quantile, prob = .025)
      W_u.CI <- apply(object$boot$weights,2,stats::quantile, prob = .975)
      
    }else{
      M_l.CI <- NULL
      M_u.CI <- NULL
      W_l.CI <- NULL
      W_u.CI <- NULL
    }
    
    TAB <- cbind('Mean '  = object$fit$mu,
                 'l.CI'   = M_l.CI,
                 'u.CI '  = M_u.CI,
                 'Weight' = object$fit$weights,
                 'l.CI'   = W_l.CI,
                 'u.CI'   = W_u.CI)
    rownames(TAB) <- as.character(1:length(object$fit$mu))
    
  }
  
  res <- list(call         = object$call,
              coefficients = TAB,
              model        = model,
              converged    = object$converged,
              round.coef   = round.coef)
  class(res) <- "summary.zcurve"
  return(res)
}

#' Prints summary object for z-curve method
#' @param x Summary of a z-curve object
#' @param ... Additional arguments
#' @method print.summary zcurve
#' @export print.summary.zcurve
#' @rawNamespace S3method(print, summary.zcurve)
#' @seealso [zcurve()]
print.summary.zcurve <- function(x, ...){
  
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat(paste(c("model: ",x$model$model, " via ", x$model$method, "\n"), collapse = ""))
  cat("\n")
  
  #stats::printCoefmat(x$coefficients, digits = 2,
  #                    cs.ind  = c(1:ncol(x$coefficients)), tst.ind = integer(), zap.ind = integer())
  
  temp_to_int  <- !rownames(x$coefficients) %in% c("Expected N", "Missing N")
  temp_coef    <- x$coefficients
  
  if(length(temp_to_int) != 0){
    temp_coef[temp_to_int,]  <- apply(as.data.frame(x$coefficients[temp_to_int,]), 2, function(p).rXd(p, X = x$round.coef))
    temp_coef[!temp_to_int,] <- round(x$coefficients[!temp_to_int,])
  }

  print(temp_coef,quote = FALSE, right = T)
  
  if(length(x$model$fit_index) > 1){
    fit_index_CI <- paste(c(", 95% CI[", .r2d(x$model$fit_index[2]), ", ", .r2d(x$model$fit_index[3]),"]"), collapse = "")
  }else{
    fit_index_CI <- NULL
  }
  cat("\n")
  
  if(x$converged){
    cat(paste(c("Model converged in ", x$model$iter, " iterations", "\n"), collapse = ""))
  }else{
    cat(paste(c("\033[0;31m", "Model did not converge in ", x$model$iter, " iterations", "\033[0m", "\n"), collapse = ""))
  }
  
  obs_proportion <- stats::prop.test(x$model$N_sig, x$model$N_all)
  cat(paste0("Fitted using ", x$model$N_used, " ", paste(x$model$input_type, collapse = " and "), "-values. ", x$model$N_all, " supplied, ", x$model$N_sig, " significant (ODR = ",  .r2d(obs_proportion$estimate), ", 95% CI [", .r2d(obs_proportion$conf.int[1]), ", ", .r2d(obs_proportion$conf.int[2]), "]).\n"))

  cat(paste(c(x$model$fit_stat," = " , .r2d(x$model$fit_index[1]), fit_index_CI, "\n"), collapse = ""))
  
}

#' Plot fitted z-curve object
#'
#' @param x Fitted z-curve object
#' @param annotation Add annotation to the plot. Defaults 
#' to \code{FALSE}.
#' @param CI Plot confidence intervals for the estimated z-curve. Defaults 
#' to \code{FALSE}.
#' @param extrapolate Scale the chart to the extrapolated area. Defaults 
#' to \code{FALSE}.
#' @param y.anno A vector of length 8 specifying the y-positions
#'   of the individual annotation lines relative to the figure's height.
#'   Defaults to \code{c(.95, .88, .78, .71, .61, .53, .43, .35)}
#' @param x.anno A number specifying the x-position of the block
#'   of annotations relative to the figure's width.
#' @param cex.anno A number specifying the size of the annotation text.
#' @param ... Additional arguments including \code{main}, \code{xlab},
#' \code{ylab}, \code{cex.axis}, \code{cex.lab}
#'
#' @method plot zcurve
#' @export plot.zcurve
#' @rawNamespace S3method(plot, zcurve)
#' 
#' @examples 
#' # simulate some z-statistics and fit a z-curve
#' z <- abs(rnorm(300,3))
#' m.EM <- zcurve(z, method = "EM", bootstrap = 100)
#' 
#' # plot the z-curve
#' plot(m.EM)
#' 
#' # add annotation text and model fit CI
#' plot(m.EM, annotation = TRUE, CI = TRUE)
#' 
#' # change the location of the annotation to the left
#' plot(m.EM, annotation = TRUE, CI = TRUE, x_text = 0)
#' @seealso [zcurve()]
plot.zcurve          <- function(x, annotation = FALSE, CI = FALSE, extrapolate = FALSE,
                                 y.anno = c(.95, .88, .78, .71, .61, .53, .43, .35), x.anno = .6, cex.anno = 1, ...){
 
  if(is.null(x$boot))CI <- FALSE
  
  additional <- list(...)
  if(is.null(additional$main)){
    main <- paste(c("z-curve (", ifelse(is.null(x$control$model), "custom", x$control$model), " via ", x$method, ")"), collapse = "")
  }else{
    main <- additional$main
  }
  if(is.null(additional$xlab)){
    xlab <- "z-scores"
  }else{
    xlab <- additional$xlab
  }
  if(is.null(additional$ylab)){
    ylab <- "Density"
  }else{
    ylab <- additional$ylab
  }
  if(is.null(additional$cex.axis)){
    cex.axis <- 1
  }else{
    cex.axis <- additional$cex.axis
  }
  if(is.null(additional$cex.lab)){
    cex.lab <- 1
  }else{
    cex.lab <- additional$cex.lab
  }
  
  # set breaks for the histogram
  br1 <- seq(x$control$a, x$control$b, .20)
  br2 <- seq(0, x$control$a, .20)
  # change the last breake to the cutpoints
  br1[length(br1)] <- x$control$b
  br2[length(br2)] <- x$control$a
  
  # get histograms
  h1 <- graphics::hist(x$data[x$data > x$control$a & x$data < x$control$b], breaks = br1, plot = F) 
  if(length(x$data[x$data < x$control$a])){
    h2 <- graphics::hist(x$data[x$data < x$control$a], breaks = br2, plot = F)
    # scale the density of nonsignificant z-scores appropriately to the first one
    h2$density <- h2$density * (x$control$a/(x$control$b - x$control$a))
    h2$density <- h2$density/(
      (length(x$data[x$data > x$control$a & x$data < x$control$b])/(x$control$b - x$control$a))
      /
        (length(x$data[x$data < x$control$a])/(x$control$a))
    )    
  }else{
    h2 <- NULL
  }

  # compute fitted z-curve density
  x_seq <- seq(0, x$control$b, .01)
  y_den <- sapply(1:length(x$fit$mu), function(i){
    x$fit$weights[i]*exp(.zdist_lpdf(x_seq, x$fit$mu[i], 1, x$control$a, x$control$b))
  })
  y_den <- apply(y_den, 1, sum)
  # and the piecewise confidence intervals
  if(CI & !is.null(x$boot)){
    y_den_boot <- sapply(1:nrow(x$boot$mu),function(b){
      y_den <- sapply(1:length(x$boot$mu[b,]), function(i){
        x$boot$weights[b,i]*exp(.zdist_lpdf(x_seq, x$boot$mu[b,i], 1, x$control$a, x$control$b))
      })
      y_den <- apply(y_den, 1, sum)
    })
    y_den_l.CI <- apply(y_den_boot, 1, stats::quantile, prob = .025)
    y_den_u.CI <- apply(y_den_boot, 1, stats::quantile, prob = .975)
  }
  
  ### setting of the axis and values for text allignment
  x_max  <- x$control$b
  x.anno <- x_max*x.anno
  
  if(extrapolate){
    y_max  <- max(c(y_den, h1$density, h2$density))
  }else{
    y_max <-  max(c(h1$density, h2$density))
  }
  
  # adjusting the height of the chart so the text is higher than the highest ploted thing in the x-range of the text
  if(annotation & CI){
    y_max <- ifelse(max(c(y_den_u.CI[x.anno < x_seq], h1$density[x.anno < h1$breaks[-length(h1$density)]])) > y_max*(y.anno[length(y.anno)] - .025),
                    max(c(y_den_u.CI[x.anno < x_seq], h1$density[x.anno < h1$breaks[-length(h1$density)]]))/(y.anno[length(y.anno)] - .025), y_max)
  }else if(annotation & !CI){
    y_max <- ifelse(max(c(y_den[x.anno < x_seq], h1$density[x.anno < h1$breaks[-length(h1$density)]])) > y_max*(y.anno[length(y.anno)] - .025),
                    max(c(y_den[x.anno < x_seq], h1$density[x.anno < h1$breaks[-length(h1$density)]]))/(y.anno[length(y.anno)] - .025), y_max)
  }
  
  
  # plot z-scores used for fitting
  graphics::plot(h1,
                 freq = FALSE, density = 0, angle = 0, border = "blue",
                 xlim = c(0, x_max),
                 ylim = c(0, y_max),
                 ylab = ylab,
                 xlab = xlab,
                 main = main,
                 cex.lab = cex.lab,
                 cex.axis = cex.axis,
                 lwd = 1, las = 1)
  # and un-used z-scores
  if(!is.null(h2)){
    graphics::par(new=TRUE)
    graphics::plot(h2,
                   freq = FALSE, density = 0, angle = 0, border ="grey30",
                   xlim = c(0, x_max),
                   ylim = c(0, y_max),
                   axes = FALSE, ann = FALSE, lwd = 1, las = 1)  
  }
  # add the density estimate if the model was estimated by density
  if(x$method == "density"){
    graphics::lines(x$fit$density$x, x$fit$density$y, lty = 1, col = "grey60", lwd = 4)
  }
  # significance line
  if(x.anno*x_max < x$control$a){
    graphics::lines(rep(x$control$a,2),                                             c(0, (min(y.anno) - .025)*y_max), col = "blue", lty = 2, lwd = 1)
    graphics::lines(rep(stats::qnorm(x$control$sig_level/2, lower.tail = FALSE),2), c(0, (min(y.anno) - .025)*y_max), col = "red",  lty = 1, lwd = 2)    
  }else{
    graphics::abline(v = x$control$a,                                             col = "blue", lty = 2, lwd = 1)
    graphics::abline(v = stats::qnorm(x$control$sig_level/2, lower.tail = FALSE), col = "red",  lty = 1, lwd = 2) 
  }
  # predicted densities
  graphics::lines(x_seq, y_den, lty = 1, col = "blue", lwd = 5)
  if(CI & !is.null(x$boot)){
    graphics::lines(x_seq, y_den_l.CI, lty = 3, col = "blue", lwd = 3)
    graphics::lines(x_seq, y_den_u.CI, lty = 3, col = "blue", lwd = 3)
  }
  # add annotation
  if(annotation){
    x_summary <- summary(x)
    
    graphics::text(x.anno, y_max*y.anno[1] , paste0("Range: ",.r2d(min(x$data))," to ",.r2d(max(x$data))),
                   adj = c(0, 0), cex = cex.anno)
    
    graphics::text(x.anno, y_max*y.anno[2] , paste0(length(x$data), " tests, ", sum(x$data >= x$control$a), " significant"),
                   adj = c(0, 0), cex = cex.anno)
    
    obs_proportion <- stats::prop.test(sum(x$data >= x$control$a), length(x$data))
    graphics::text(x.anno, y_max*y.anno[3] , paste0("Observed discovery rate:"),
                   adj = c(0, 0), cex = cex.anno)
    graphics::text(x.anno, y_max*y.anno[4] , paste0(.r2d(obs_proportion$estimate), "  95% CI [", .r2d(obs_proportion$conf.int[1]), " ,",
                                                    .r2d(obs_proportion$conf.int[2]), "]"),
                   adj = c(0, 0), cex = cex.anno)
    
    if(!is.null(x$boot)){
      graphics::text(x.anno, y_max*y.anno[5] , paste0("Expected discovery rate:"),
                     adj = c(0, 0), cex = cex.anno)
      graphics::text(x.anno, y_max*y.anno[6] , paste0(.r2d(x_summary$coefficients["EDR","Estimate"]), "  95% CI [", .r2d(x_summary$coefficients["EDR","l.CI"]), " ,",
                                                      .r2d(x_summary$coefficients["EDR","u.CI"]), "]"),
                     adj = c(0, 0), cex = cex.anno)
      
      graphics::text(x.anno, y_max*y.anno[7] , paste0("Expected replicability rate:"),
                     adj = c(0, 0), cex = cex.anno)
      graphics::text(x.anno, y_max*y.anno[8] , paste0(.r2d(x_summary$coefficients["ERR","Estimate"]), "  95% CI [", .r2d(x_summary$coefficients["ERR","l.CI"]), " ,",
                                                      .r2d(x_summary$coefficients["ERR","u.CI"]), "]"),
                     adj = c(0, 0), cex = cex.anno)
    }else{
      graphics::text(x.anno, y_max*y.anno[5] , paste0("Expected discovery rate:"),
                     adj = c(0, 0), cex = cex.anno)
      graphics::text(x.anno, y_max*y.anno[6] , paste0(.r2d(x_summary$coefficients["EDR","Estimate"])),
                     adj = c(0, 0), cex = cex.anno)
      
      graphics::text(x.anno, y_max*y.anno[7] , paste0("Expected replicability rate:"),
                     adj = c(0, 0), cex = cex.anno)
      graphics::text(x.anno, y_max*y.anno[8] , paste0(.r2d(x_summary$coefficients["ERR","Estimate"])),
                     adj = c(0, 0), cex = cex.anno)
    }

    
  }
  
}

#' @title Z-scores from subset of original studies featured in OSC 2015 
#' reproducibility project
#' 
#' @description The dataset contains z-scores from subset of original
#'  studies featured in psychology reproducibility project 
#'  \insertCite{osc}{zcurve}. Only z-scores from studies with unambiguous 
#' original outcomes are supplied (eliminating 7 studies with marginally 
#' significant results). The real replication rate for those studies is 
#' 35/90 (the whole project reports 36/97).
#'
#' @format A vector with 90 observations
#' 
#' @references
#' \insertAllCited{}
"OSC.z"

# cleaning
.onUnload <- function (libpath) {
  library.dynam.unload("zcurve", libpath)
}
