# wrapper
.zcurve_density          <- function(z_sig, control) {
  
  if(control$version == 2){
    temp_fit <- .zcurve_density_ver2(z_sig, control)
  }else if(control$version == 1){
    temp_fit <- .zcurve_density_ver1(z_sig, control)
  }
  
  return(temp_fit)
  
}
.zcurve_density_boot     <- function(z_sig, control, bootstrap){
  
  temp_ncol <- ifelse(control$version == 2, length(control$mu), control$K)
  results <-  list(
    "mu"        = matrix(NA, ncol = temp_ncol, nrow = bootstrap),
    "weights"   = matrix(NA, ncol = temp_ncol, nrow = bootstrap),
    "N_fit"     = rep(NA, times = bootstrap),
    "objective" = rep(NA, times = bootstrap),
    "iter"      = rep(NA, times = bootstrap),
    "FDR_max"   = rep(NA, times = bootstrap)
    # probably not useful to save individual densities from the bootstrap
    #"density"   = list(
    #  "x" = NULL,
    #  "y" = NULL
    #)
  )
  
  for(i in 1:bootstrap){
    z_sig_boot <- sample(z_sig, replace = TRUE)
    
    if(control$version == 2){
      temp_fit <- .zcurve_density_ver2(z_sig_boot, control)
      results$FDR_max[i]   <- temp_fit$FDR_max
    }else if(control$version == 1){
      temp_fit <- .zcurve_density_ver1(z_sig_boot, control)
    }
    
    
    results$mu[i,]       <- temp_fit$mu
    results$weights[i,]  <- temp_fit$weights
    results$N_fit[i]     <- temp_fit$N_fit
    results$objective[i] <- temp_fit$objective
    results$iter[i]      <- temp_fit$iter

  }
  
  return(results)
}



# original z-curve1.0
.zcurve_density_ver1         <- function(z_sig, control) {
  
  ncomp  <- control$K
  bw <- control$bw
  cv <- control$a
  Z  <- z_sig

  augZ = c(subset(Z,Z<Inf),2*cv-subset(Z,Z<Inf)) 	#5. Augmented Z for density estimation to avoid asymtote to zero at cv.
  DensityEstimate = stats::density(augZ,n=100,bw=bw,from=1.96,to=6) #6. Get Densities using Kernel.Density.Function
  dens.Z = DensityEstimate$x;  #7. x-axis values of Density Slices (z-scores)
  dens.obs = DensityEstimate$y; #8. Observed Densities for the Density Slices
  dens.obs = dens.obs/(sum(dens.obs)) #9. Sum of Densities equals 1  
  slices = length(dens.Z) #10. Number of "Slices" of the Density Distribution
  slices.width = dens.Z[2] - dens.Z[1] #11. Width of a "Slice" of the Density Distribution
  ### This is the actual estimation function that fits estimated density to observed density
  
  ### Prepare start values and limits for nlminb package that fits z-curve to the data
  startval = c(1,2,3,1/3,1/3,1/3) #25 Starting Values for Means (1,2,3) Starting Values for weights (1/3)
  lowlim = rep(0,6) #26 lower limit for Means and Weights, all 0
  highlim = c(6,6,6,1,1,1) #27 upper limit for Means = 6, upper limit for weights = 1
  ### Execute the fit.zcurve function to get estimates
  auto = suppressWarnings(stats::nlminb(startval,.zcurve_density_ver1_fit,lower=lowlim,upper=highlim,
                                        ncomp = ncomp, dens.Z = dens.Z, dens.obs = dens.obs,
                                        control=list(eval.max=control$max_eval,
                                                     iter.max = control$max_iter,
                                                     rel.tol = control$criterion)
                                        )
                          )
  #28 nlminb searches for Means and Weights that minimize the fit criterion
  ### Get the Estimated Means and Weights
  Z.Means = auto$par[1:ncomp] #29 get the final Means
  Z.w = auto$par[(ncomp+1):(2*ncomp)] #30 Get the final Weights
  Z.w = Z.w/sum(Z.w) #31 Weights are positive and sum to one.
  
  mean    = Z.Means
  weights = Z.w
  
  # scaling for plotting 
  bar.width   = DensityEstimate$x[2] - DensityEstimate$x[1]
  Z.Density.Y = DensityEstimate$y/(sum(DensityEstimate$y*bar.width))
  
  return(
    list(
      "mu"        = Z.Means,
      "weights"   = Z.w,
      "N_fit"     = sum(z_sig < control$b),
      "objective" = auto$objective,
      "converged" = auto$convergence,
      "message"   = auto$message,
      "iter"      = auto$iterations,
      "density"   = list(
        "x" = DensityEstimate$x,
        "y" = Z.Density.Y
      )
    )
  )
  
}
.zcurve_density_ver1_fit     <- function(parameter, ncomp, dens.Z, dens.obs) {   #12 repeated until best fit is reached 
  Z.Means = parameter[1:ncomp] #13 These are the estimated means while approximating observed density
  Z.weights = abs(parameter[(ncomp+1):(2*ncomp)]) #14 These are the estimated weights 
  Z.weights = Z.weights/sum(Z.weights) #15 Weights are scaled to add to 1
  
  dens = c()  #16 variable that stores the estimated densities for each component
  for(j in 1:ncomp) { #17 Do for each component 
    dcomp = stats::dnorm(dens.Z-Z.Means[j]) #18 get the densities for the z-scores of the observed density distribution
    dcomp = dcomp/sum(dcomp)  #19 scale them to add up to 1
    dcomp = dcomp*Z.weights[j] #20 weight them according to the estimated weight
    dens = rbind(dens,dcomp) ##21 store results in a matrix 
  }  # End of For loop 
  dens.est = colSums(dens) #22 get the estimated density as the sum of the weighted densities of the 3 components
  MeanAbsError = mean(abs(dens.est - dens.obs)) #23 Mean absolute difference is the fit criterion (smaller values = better fit)
  return(MeanAbsError) #24 return fit value to the fitting function
}

#' @name control_density_v1
#' @title Control settings for the original z-curve density algorithm
#' @description All settings are passed to the density fitting
#' algorithm. All unspecified settings are set to the default value.
#' Setting \code{model = "KD1"} sets all settings to the default
#' value irrespective of any other setting and fits z-curve as described 
#' in \insertCite{zcurve1;textual}{zcurve}.
#' 
#' @param version Set to \code{1} to fit the original version of z-curve. 
#' Defaults to \code{2} = the updated version of z-curve. For its settings 
#' page go to [control_density].
#' @param model A type of model to be fitted, defaults to \code{"KD1"}
#' (the only possibility)
#' @param sig_level An alpha level of the test statistics, defaults to
#' \code{.05}
#' @param a A beginning of fitting interval, defaults to
#' \code{qnorm(sig_level/2,lower.tail = F)}
#' @param b An end of fitting interval, defaults to \code{6}
#' @param K Number of mixture components, defaults to \code{3}
#' @param max_iter A maximum number of iterations for the \link[stats]{nlminb} 
#' optimization for fitting mixture model, defaults to \code{150}
#' @param max_eval A maximum number of evaluation for the \link[stats]{nlminb} 
#' optimization for fitting mixture model, defaults to \code{300}
#' @param criterion A criterion to terminate \link[stats]{nlminb} optimization,
#' defaults to \code{1e-10}
#' @param bw A bandwidth of the kernel density estimation, defaults to \code{"nrd0"}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @examples # to increase the number of iterations
#' ctrl <- list(
#'    version   = 1,
#'    max_iter  = 300
#' )
#' \dontrun{zcurve(z, algorithm = "density", control = ctrl)}
#' 
#' @seealso [zcurve()], [control_density], [control_EM]
NULL


# z-curve 2.0 (KD2)
.zcurve_density_ver2      <- function(z_sig, control) {

  z.val.input <- z_sig
  z.val.input[z.val.input > control$MAX.INP.Z] = control$MAX.INP.Z

  Z.INT = z.val.input[z.val.input >= control$a & z.val.input <= control$b + 1]

  z.extreme = sum(z.val.input > control$b)/sum(z.val.input > control$a)

  densy = .zcurve_density_get_densities(Z.INT, z.val.input, control)

  Z.Density.X <- densy[,1]
  Z.Density.Y <- densy[,2]

  control$SLOPE = stats::coef(stats::lm(Z.Density.Y[Z.Density.X < control$a+1]  ~ Z.Density.X[Z.Density.X < control$a+1]))[1]

  n.bars = length(Z.Density.X)

  bar.width = Z.Density.X[2] - Z.Density.X[1]

  ### get the densities for each interval and each non-centrality parameter
  Dens	= c()
  for(i in 1:n.bars) {
    for (j in 1:length(control$mu)) {
      Dens = c(Dens,.zdist_pdf(Z.Density.X[i],control$mu[j],control$sigma,control$a,control$b))
    }
  }
  Dens = matrix(Dens,length(control$mu),byrow=FALSE)

  Dens = Dens/(rowSums(Dens) * bar.width)

  para.val = .zcurve_density_get_weights_free(control, Dens = Dens, n.bars = n.bars,
                                              Z.Density.Y = Z.Density.Y, Z.Density.X = Z.Density.X)
  # control$SLOPE = para.val[1] probably redundant
  WZ0 = para.val$weights[1]
  fit.free = para.val$objective
  FDR.RES = c(WZ0,NA)
  FDR.RES = FDR.RES*(1-z.extreme)

  #para.est = Compute.Power(c(para.val$mu, para.val$weights),z.extreme)
  #para.est

  W = para.val$weights
  #W
  # precision = .2
  if (control$compute_FDR) FDR.RES[2] = .zcurve_density_get_weights_fixed(z.val.input    = z.val.input,
                                                                          W              = W,
                                                                          fit.free       = fit.free,
                                                                          precision      = .2,
                                                                          n.bars         = n.bars,
                                                                          Z.Density.X    = Z.Density.X,
                                                                          Z.Density.Y    = Z.Density.Y,
                                                                          Dens           = Dens,
                                                                          control        = control)
  #FDR.RES[2]
  # precision = control$precision_FDR
  W[1] = max(FDR.RES)
  if (control$compute_FDR) FDR.RES[2] = .zcurve_density_get_weights_fixed(z.val.input    = z.val.input,
                                                                          W              = W,
                                                                          fit.free       = fit.free,
                                                                          precision      = control$precision_FDR,
                                                                          n.bars         = n.bars,
                                                                          Z.Density.X    = Z.Density.X,
                                                                          Z.Density.Y    = Z.Density.Y,
                                                                          Dens           = Dens,
                                                                          control        = control)
  #FDR.RES

  #res = c(para.est,FDR.RES)
  #names(res) = c("ERR","EDR","Weight0","Max.FDR")
  #res

  return(
    list(
      "mu"        = para.val$mu,
      "weights"   = para.val$weights,
      "N_fit"     = sum(z_sig < control$b),
      "objective" = para.val$objective,
      "converged" = para.val$converged,
      "message"   = para.val$message,
      "iter"      = para.val$iter,
      "FDR_max"   = FDR.RES[2],
      "density"   = list(
        "x" = Z.Density.X,
        "y" = Z.Density.Y
      )
    )
  )

}


#' @name control_density
#' @title Control settings for the z-curve 2.0 density algorithm
#' @description All settings are passed to the density fitting
#' algorithm. All unspecified settings are set to the default value.
#' Setting \code{model = "KD2"} sets all settings to the default
#' value irrespective of any other setting and fits z-curve as 
#' describe in \insertCite{zcurve2;textual}{zcurve}. In order to fit the 
#' z-curve 1.0 density algorithm, set \code{model = "KD1"} and go to 
#' [control_density_v1]
#' 
#' @param version Which version of z-curve should be fitted. Defaults to
#'  \code{2} = z-curve 2.0. Set to \code{1} in order to fit the original 
#'  version of z-curve. For its settings page go to [control_density_v1].
#' @param model A type of model to be fitted, defaults to \code{"KD2"}
#' (another possibility is \code{"KD1"} for the original z-curve 1.0, see
#' [control_density_v1] for its settings)
#' @param sig_level An alpha level of the test statistics, defaults to
#' \code{.05}
#' @param a A beginning of fitting interval, defaults to
#' \code{qnorm(sig_level/2,lower.tail = F)}
#' @param b An end of fitting interval, defaults to \code{6}
#' @param mu Means of the components, defaults to \code{seq(0,6,1)}
#' @param sigma A standard deviation of the components, "Don't touch this"
#' \- Ulrich Schimmack, defaults to \code{1}
#' @param theta_min Lower limits for weights, defaults to
#' \code{rep(0,length(mu))}
#' @param theta_max Upper limits for weights, defaults to
#' \code{rep(1,length(mu))}
#' @param max_iter A maximum number of iterations for the \link[stats]{nlminb} 
#' optimization for fitting mixture model, defaults to \code{150}
#' @param max_eval A maximum number of evaluation for the \link[stats]{nlminb} 
#' optimization for fitting mixture model, defaults to \code{1000}
#' @param criterion A criterion to terminate \link[stats]{nlminb} optimization,
#' defaults to \code{1e-03}
#' @param bw A bandwidth of the kernel density estimation, defaults to \code{.10}
#' @param aug Augment truncated kernel density, defaults to \code{TRUE}
#' @param aug.bw A bandwidth of the augmentation, defaults to \code{.20}
#' @param n.bars A resolution of density function, defaults to \code{512}
#' @param density_dbc Use \link[evmix]{bckden} to estimate a truncated kernel density,
#' defaults to \code{FALSE}, in which case \link[stats]{density} is used
#' @param compute_FDR Whether to compute FDR, leads to noticeable increase in
#' computation, defaults to \code{FALSE}
#' @param criterion_FDR A criterion for estimating the maximum FDR, defaults
#' to \code{.02}
#' @param criterion_FDR_dbc A criterion for estimating the maximum FDR using
#' the \link[evmix]{bckden} function, defaults to \code{.01}
#' @param precision_FDR A maximum FDR precision, defaults to \code{.05}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @examples # to decrease the criterion and increase the number of iterations
#' ctrl <- list(
#'    max_iter  = 300,
#'    criterion = 1e-4
#' )
#' \dontrun{zcurve(z, algorithm = "density", control = ctrl)}
#' 
#' @seealso [zcurve()], [control_density_v1], [control_EM]
NULL

.zcurve_density.control  <- function(control){
  if(is.null(control)){
    control$version         <- 2
    control$sig_level       <- .05
    control$sig_level_Z     <- stats::qnorm(control$sig_level/2,lower.tail = F)
    control$a               <- stats::qnorm(control$sig_level/2,lower.tail = F)  	# Beginning of fitting interval
    control$b               <- 6  	          	# End of fitting interval
    control$mu              <- seq(0,6,1)  	    # means of the components
    control$sigma           <- 1                # Don't touch this!!! # Change Standard Deviation of Normals
    control$theta_min       <- rep(0,length(control$mu))	# set lower limit for weights,  default = 0
    control$theta_max       <- rep(1,length(control$mu))	# set upper limits for weights, default = 1
    control$max_iter        <- 150              # settings for the nlminb agortihm for fitting mixture model
    control$max_eval        <- 1000             # settings for the nlminb agortihm for fitting mixture model
    control$criterion       <- 1e-03           	# Criterion to terminate nlminb
    control$bw              <- .10              # Bandwidth of Kernal Density
    control$aug             <- TRUE		          # Augment truncated Kernal Density
    control$aug.bw          <- .20  		        # Augment Bandwidth
    control$n.bars          <- 512   		        # resolution of density function (doesn't seem to matter much)
    control$density_dbc     <- FALSE			      # USE dbckden function to truncate Kernal Density
    control$criterion_FDR     <- .02		        # criterion for maximum FDR
    control$criterion_FDR_dbc <- .01	          # criterion for maximum FDR using density_dbc function
    control$precision_FDR   <- .05              # Maximum FDR precision (low precision slows things down)
    control$compute_FDR     <- FALSE	          # Compute Maximum FDR, Slows Things Down Considerably

    control$model           <- "KD2"
    
    # check whether this is needed
    control$MAX.INP.Z       <- 100	            # values greater than MAX.INP.Z will be set to maximum

    ### probably to be removed
    control$PLOT = F
    control$FDR.PLOT = F
    #"FDR.PLOT" = FALSE	   # Make a screen plot of the FDR
    #"PLOT" = FALSE		     # Show Fitting of Density Distribution

    ### unused in the code
    #"SLOPE.crit" = 1		# Slope Criterion to set EDR estimate to NA
    #"USE.SLOPE" = FALSE	# Use Slope Criterion to Exclude EDR estimates

    ### unused and should be elsewhere
    #"BOOT.FDR" = TRUE		# Do a Bootstrap for the Max. FDR, really slow
    #"BOOT" = FALSE		# Bootstrap or No.Bootstrap Computation of Power
    #"boot.iter" = 0		#  How many bootstraps for CI; 0 = no CI

    return(control)
  }
  if(!is.null(control$model)){
    if(control$model == "KD2"){
      control$version         <- 2
      control$sig_level       <- .05
      control$sig_level_Z     <- stats::qnorm(control$sig_level/2,lower.tail = F)
      control$a               <- stats::qnorm(control$sig_level/2,lower.tail = F)  	# Beginning of fitting interval
      control$b               <- 6  	          	# End of fitting interval
      control$mu              <- seq(0,6,1)  	    # means of the components
      control$sigma           <- 1                # Don't touch this!!! # Change Standard Deviation of Normals
      control$theta_min       <- rep(0,length(control$mu))	# set lower limit for weights,  default = 0
      control$theta_max       <- rep(1,length(control$mu))	# set upper limits for weights, default = 1
      control$max_iter        <- 150              # settings for the nlminb agortihm for fitting mixture model
      control$max_eval        <- 1000             # settings for the nlminb agortihm for fitting mixture model
      control$criterion       <- 1e-03           	# Criterion to terminate nlminb
      control$bw              <- .10              # Bandwidth of Kernal Density
      control$aug             <- TRUE		          # Augment truncated Kernal Density
      control$aug.bw          <- .20  		        # Augment Bandwidth
      control$n.bars          <- 512   		        # resolution of density function (doesn't seem to matter much)
      control$density_dbc     <- FALSE			      # USE dbckden function to truncate Kernal Density
      control$criterion_FDR     <- .02		        # criterion for maximum FDR
      control$criterion_FDR_dbc <- .01	          # criterion for maximum FDR using density_dbc function
      control$precision_FDR   <- .05              # Maximum FDR precision (low precision slows things down)
      control$compute_FDR     <- FALSE	          # Compute Maximum FDR, Slows Things Down Considerably
      
      control$model           <- "KD2"
      
      # check whether this is needed
      control$MAX.INP.Z       <- 100	            # values greater than MAX.INP.Z will be set to maximum
      
      ### probably to be removed
      control$PLOT = F
      control$FDR.PLOT = F
      
      return(control)
    }
    if(control$model == "KD1"){
      control$version         <- 1
      control$sig_level       <- .05
      control$sig_level_Z     <- stats::qnorm(control$sig_level/2,lower.tail = F)
      control$a               <- stats::qnorm(control$sig_level/2,lower.tail = F)  	# Beginning of fitting interval
      control$b               <- 6  	          	# End of fitting interval
      control$K               <- 3
      control$max_iter        <- 150              # settings for the nlminb agortihm for fitting mixture model
      control$max_eval        <- 300              # settings for the nlminb agortihm for fitting mixture model
      control$criterion       <- 1e-10           	# Criterion to terminate nlminb
      control$bw              <- "nrd0"           # Bandwidth of Kernal Density
      control$model           <- "KD1"
      
      return(control)
    }
  }
  if(!is.null(control$version)){
    if(control$version == 1){
      if(is.null(control$sig_level)){
        control$sig_level       <- .05
      }
      if(is.null(control$sig_level_Z)){
        control$sig_level_Z     <- stats::qnorm(control$sig_level/2,lower.tail = F)
      }
      if(is.null(control$a)){
        control$a               <- stats::qnorm(control$sig_level/2,lower.tail = F)  	# Beginning of fitting interval
      }
      if(is.null(control$b)){
        control$b               <- 6  	          	# End of fitting interval
      }
      if(is.null(control$K)){
        control$K               <- 3
      }
      if(is.null(control$max_iter)){
        control$max_iter        <- 150              # settings for the nlminb agortihm for fitting mixture model
      }
      if(is.null(control$max_eval)){
        control$max_eval        <- 300             # settings for the nlminb agortihm for fitting mixture model
      }
      if(is.null(control$criterion)){
        control$criterion       <- 1e-10           	# Criterion to terminate nlminb
      }
      if(is.null(control$bw)){
        control$bw              <- "nrd0"           # Bandwidth of Kernal Density
      }
      return(control)
    }
  }
  # individual parameter settings
  if(is.null(control$version)){
    control$version         <- 2
  }
  if(is.null(control$sig_level)){
    control$sig_level       <- .05
  }
  if(is.null(control$sig_level_Z)){
    control$sig_level_Z     <- stats::qnorm(control$sig_level/2,lower.tail = F)
  }
  if(is.null(control$a)){
    control$a               <- stats::qnorm(control$sig_level/2,lower.tail = F)  	# Beginning of fitting interval
  }
  if(is.null(control$b)){
    control$b               <- 6  	          	# End of fitting interval
  }
  if(is.null(control$mu)){
    control$mu              <- seq(0,6,1)  	    # means of the components
  }
  if(is.null(control$sigma)){
    control$sigma           <- 1                # Don't touch this!!! # Change Standard Deviation of Normals
  }
  if(is.null(control$theta_min)){
    control$theta_min       <- rep(0,length(control$mu))	# set lower limit for weights,  default = 0
  }
  if(is.null(control$theta_max)){
    control$theta_max       <- rep(1,length(control$mu))	# set upper limits for weights, default = 1
  }
  if(is.null(control$max_iter)){
    control$max_iter        <- 150              # settings for the nlminb agortihm for fitting mixture model
  }
  if(is.null(control$max_eval)){
    control$max_eval        <- 1000             # settings for the nlminb agortihm for fitting mixture model
  }
  if(is.null(control$criterion)){
    control$criterion       <- 1e-03           	# Criterion to terminate nlminb
  }
  if(is.null(control$bw)){
    control$bw              <- .10              # Bandwidth of Kernal Density
  }
  if(is.null(control$aug)){
    control$aug             <- TRUE		          # Augment truncated Kernal Density
  }
  if(is.null(control$aug.bw)){
    control$aug.bw          <- .20  		        # augation Bandwidth
  }
  if(is.null(control$n.bars)){
    control$n.bars          <- 512   		        # resolution of density function (doesn't seem to matter much)
  }
  if(is.null(control$density_dbc)){
    control$density_dbc     <- FALSE			      # USE dbckden function to truncate Kernal Density
  }
  if(is.null(control$criterion_FDR)){
    control$criterion_FDR     <- .02		        # criterion for maximum FDR
  }
  if(is.null(control$criterion_FDR_dbc)){
    control$criterion_FDR_dbc <- .01	          # criterion for maximum FDR using density_dbc function
  }
  if(is.null(control$precision_FDR)){
    control$precision_FDR   <- .05              # Maximum FDR precision (low precision slows things down)
  }
  if(is.null(control$compute_FDR)){
    control$compute_FDR     <- FALSE	          # Compute Maximum FDR, Slows Things Down Considerably
  }
  if(is.null(control$model)){
    control$model           <- NULL
  }
  # check whether this is needed
  if(is.null(control$MAX.INP.Z)){
    control$MAX.INP.Z       <- 100	            # values greater than MAX.INP.Z will be set to maximum
  }
  ### probably to be removed
  if(is.null(control$PLOT)){
    control$PLOT = F
  }
  if(is.null(control$FDR.PLOT)){
    control$FDR.PLOT = F
  }
  return(control)
}

#### additional functions ####
# not used anymore
#########################################################################
### This Function Computes Power from Weights and Non-Centrality Parameters
#########################################################################
#Compute.Power = function(para.val,z.extreme) {
#
#  ### the input weights based on z.curve method
#  ### these are the weights based on the a value
#  ### a could be 1.96 (all significant)
#  ### but it can also be other values
#  ### Therefore the weights cannot be directly used
#  ### to estimate power/replicability
#  w.inp = para.val[(length(control$mu)+1):(length(control$mu)*2)]#
#
#  pow = pnorm(control$mu,control$sig_level_Z) + pnorm(-control$mu,control$sig_level_Z)#
#
#  ### this gives the power with the a z-score as the criterion value
#  ### this power is used as a weight to get the weights for the full distribution
#  ### using Jerry's insight that weight before selection is weight after selection divided by power
#  w = pnorm(control$mu,control$a) + pnorm(-control$mu,control$a)
#  round(w,3)#
#
#  ### now we compute the weights before selection (w.all)
#  ### once we have the weights, we devided by sum of all weights
#  ### so that they add up to 1
#  w.all = w.inp / w
#  w.all = w.all / sum(w.all)
#
#  ### now we are ready to compute the weights after selection for significance
#  ### using Jerry's fomrula in reverse going from before selection to after selection
#  ### by multiplying by power (w)
#  ### again all the weights are standardized by dividing by the sum of all weights
#
#  w.sig = w.all * pow
#  w.sig = w.sig / sum(w.sig)
#  w.sig
#
#  ### compute ERR
#  ### this is easy, replicabilty is simply the weighted sum of power
#  ### using the weights after selection for significance, w.sig
#  ERR = sum(pow*w.sig)
#
#  ### than the maximum value used (default z > 6)
#  ERR = ERR*(1 - z.extreme) + z.extreme
#
#  ### compute average power for all results, including estimated file drawer
#  ### this is also easy, here average power before selection is computed
#  ### as the weighted average of power using the weights before selection, w.all
#
#  w.ext = c(w.sig*(1-z.extreme),z.extreme)
#  pow.ext = c(pow,1)
#  EDR = 1/sum(w.ext/pow.ext)
#
#  ### res stores the results to be past back from the function
#  res = c(ERR,EDR)
#
#  return(res)
#
#}



#######################################################
### Fitting ZCurve for FDR ESTIMATE
#######################################################
.zcurve_density_get_weights_fixed = function(z.val.input, W, fit.free, precision,
                                             n.bars, Z.Density.X, Z.Density.Y, Dens, control){


  if (control$FDR.PLOT) {

    fit = c()

    W.set = seq(0,1,control$precision_FDR)

    for (WZ0 in W.set) {

      theta_min = rep(0,length(control$mu))
      theta_max = rep(1,length(control$mu))

      startval = (control$theta_min+control$theta_max)/2
      startval = startval/sum(startval)
      startval[1] = 1

      theta_min[1] = 0
      theta_max[1] = 0

      ### start the estimation process
      auto = stats::nlminb(startval,.zcurve_density_fitting_fixed,lower=theta_min,upper=theta_max,
                           WZ0 = WZ0, n.bars = n.bars,
                           Z.Density.X = Z.Density.X, Z.Density.Y = Z.Density.Y, Dens = Dens,
                           PLOT = control$PLOT)
      fit = c(fit,auto$objective)

    }

    graphics::plot((W.set-1)/10,fit,ylim=c(0,max(fit)+.05),xlab="Percentage of False Positives",ylab="Root Mean Square Discrepancy of Densities", col="red",
         pch=16, cex=2)
    graphics::abline(h=fit.free,col="blue",lwd=2,lty=2)
    graphics::abline(h=fit.free+crit,lty=2,col="red")
    # windows()

    crit = fit.free + control$criterion_FDR
    crit
    if (control$density_dbc) crit = fit.free + control$criterion_FDR_dbc
    MAX.FDR = max(W.set[fit < crit])
    MAX.FDR

  } else {

    WZ0 = trunc(W[1]/precision)*precision

    crit = control$criterion_FDR + fit.free
    if (control$density_dbc) crit = control$criterion_FDR_dbc + fit.free

    fit.z0 = 0

    while (fit.z0 < crit & WZ0 <= 1) {

      theta_min = rep(0,length(control$mu))
      theta_max = rep(1,length(control$mu))

      startval = W+.10
      startval = startval/sum(startval)

      theta_min[1] = 0
      theta_max[1] = 0

      auto = stats::nlminb(startval,.zcurve_density_fitting_fixed,
                    control=list(rel.tol = control$criterion), lower=theta_min,upper=theta_max,
                    WZ0 = WZ0, n.bars = n.bars, Z.Density.Y = Z.Density.Y, PLOT = control$PLOT)
      auto$par

      W = (auto$par/sum(auto$par))*(1-WZ0)
      W[1] = WZ0
      W
      fit.z0 = auto$objective
      if (fit.z0 < crit) WZ0 = WZ0 + precision
    }
    MAX.FDR = WZ0 - precision
  }


  return(MAX.FDR)

}
.zcurve_density_fitting_fixed     = function(theta, WZ0, n.bars, Z.Density.X, Z.Density.Y, Dens, PLOT){

  ### get the weights and rescale
  theta = theta/sum(theta)*(1-WZ0)
  weight = c(WZ0,theta)

  ### compute the new estimated density distribution
  z.est = c()
  for (i in 1:n.bars) z.est[i] = sum(Dens[,i]*weight)

  ### compare to observed density distribution
  rmse = sqrt(mean((z.est-Z.Density.Y)^2))

  ### showing the fitting of the function in a plot
  if(PLOT) {

    rval = stats::runif(1)
    if (rval > .9) {
      graphics::lines(Z.Density.X,z.est,lty=1,col="red1",ylim=c(0,1),)
      graphics::points(Z.Density.X,z.est,pch=20,col="red1",ylim=c(0,1),)

    }

  }

  ### return value to optimization function
  return(rmse)
}

#######################################################
### Get Weights of Non-Central Z-scores
#######################################################
### This function fits an observed distribution of z-scores to a multi-model mixture model
.zcurve_density_get_weights_free = function(control, Dens, n.bars, Z.Density.Y, Z.Density.X){

  startval = rep(1/length(control$mu),length(control$mu))
  startval[1] = 1
  startval = startval/sum(startval)

  auto = stats::nlminb(startval,.zcurve_density_fitting_free,
                       Dens = Dens, n.bars = n.bars, Z.Density.Y = Z.Density.Y, Z.Density.X = Z.Density.X,
                       lower = control$theta_min, upper = control$theta_max, PLOT = control$PLOT,
                       control = list(eval.max=control$max_eval, iter.max = control$max_iter))

  fit.free = auto$objective

  WT = auto$par
  WT = WT/sum(WT)

  # res = c(control$SLOPE,control$mu,WT,fit.free)
  # names(res) = c("SLOPE",rep("mu",length(control$mu)),rep("Weight",length(WT)),"Fit.Free")

  res <- list(
    "slope"     = control$SLOPE,
    "mu"        = control$mu,
    "weights"   = WT,
    "objective" = auto$objective,
    "iter"      = auto$iterations,
    "converged" = auto$convergence == 0,
    "message"   = auto$message
  )

  return(res)

}
### THIS IS THE FUNCTION THAT COMPARES OBSERVED TO PREDICTED Z-VALUE DISTRIBUTIONS
.zcurve_density_fitting_free     = function(theta, Dens, n.bars, Z.Density.Y, Z.Density.X, PLOT){

  ### get the weights and rescale
  weight = theta
  weight = weight/sum(weight)

  ### compute the new estimated density distribution
  z.est = c()
  for (i in 1:n.bars) z.est[i] = sum(Dens[,i]*weight)

  ### compare to observed density distribution
  rmse = sqrt(mean((z.est-Z.Density.Y)^2))


  ### showing the fitting of the function in a plot
  if(PLOT==TRUE) {

    if (stats::runif(1) < .1) {
      graphics::plot(Z.Density.X,Z.Density.Y,type='l',ylim=c(0,1),xlab='Z')
      graphics::lines(Z.Density.X,z.est,lty=1,col="red1",ylim=c(0,1),)
      graphics::points(Z.Density.X,z.est,pch=20,col="red1",ylim=c(0,1),)
      Sys.sleep(1)
    }

  }


  ### return value to optimization function
  return(rmse)

}

##############################################
### Get Densities
##############################################
.zcurve_density_get_densities = function(Z.INT, z.val.input, control){

  ### find the maximum z-score. This is only needed if the maximum z-score is below b

  max.z = control$b
  if (max(Z.INT) < max.z) max.z = max(Z.INT)

  if (control$density_dbc) {
    Z.Density.X = seq(control$a,control$b,.01)-control$a
    xx = Z.INT-control$a
    Z.Density.Y = evmix::dbckden(Z.Density.X,xx,bw=control$bw,bcmethod="reflect")
    Z.Density.X = Z.Density.X + control$a
  } else {
    if (control$aug) {
      if (control$a >= control$sig_level_Z + 2*control$bw) {
        AUG = z.val.input[z.val.input > control$a - 2*control$bw & z.val.input < control$a]
      } else {
        AUG = c()
        n.AUG = round(length(Z.INT[Z.INT > control$a & Z.INT < control$a+control$aug.bw]))
        if (n.AUG > 0) AUG = seq(control$a-control$aug.bw,control$a-.01,control$aug.bw/n.AUG)
      }

      Z.INT.USE = c(Z.INT,AUG)

    } else {
      Z.INT.USE = Z.INT[Z.INT > control$a & Z.INT <= max.z + 1]
    }


    Z.Density = stats::density(Z.INT.USE,n=control$n.bars,bw=control$bw,from=control$a,to=max.z)
    Z.Density.X = cbind(Z.Density$x,Z.Density$y)[Z.Density$x > control$a & Z.Density$x < max.z,1]
    Z.Density.Y = cbind(Z.Density$x,Z.Density$y)[Z.Density$x > control$a & Z.Density$x < max.z,2]

  } # End of density_dbc

  bar.width = Z.Density.X[2] - Z.Density.X[1]
  Z.Density.Y = Z.Density.Y/(sum(Z.Density.Y*bar.width))

  densy = cbind(Z.Density.X,Z.Density.Y)

  return(densy)

}  ### End of Get Densities
