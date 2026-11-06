context("z-curve")
skip_on_cran()

test_that("z-curve EM can be fitted and reproduces OSC results", {

  # set seed for reproducibility 
  set.seed(666)
  
  # fit the EM method
  fit.EM  <- suppressWarnings(zcurve(OSC.z))

  # print
  expect_equal(
    capture_output_lines(print(fit.EM), print = TRUE, width = 150),
    c("Call:"               ,
      "zcurve(z = OSC.z)"   ,
      ""                    , 
      "Estimates:"          , 
      "      ERR       EDR ",
      "0.6180361 0.3905984 "    
    ))
  
  # basic summary
  expect_equal(
    capture_output_lines(summary(fit.EM), print = TRUE, width = 150),
    c("Call:"                                                                                   ,
      "zcurve(z = OSC.z)"                                                                       ,
      ""                                                                                        ,
      "model: EM via EM"                                                                        ,
      ""                                                                                        ,
      "    Estimate  l.CI  u.CI"                                                                ,
      "ERR    0.618 0.447 0.743"                                                                ,
      "EDR    0.391 0.070 0.701"                                                                ,
      ""                                                                                        ,
      "Model converged in 27 + 783 iterations"                                                  ,
      "Fitted using 73 z-values. 90 supplied, 85 significant (ODR = 0.94, 95% CI [0.87, 0.98]).",
      "Q = -60.61, 95% CI[-72.24, -46.24]"       
    ))
    
  # additional summary
  expect_equal(
    capture_output_lines(summary(fit.EM, all = TRUE), print = TRUE, width = 150),
    c("Call:"                                                                                   ,
      "zcurve(z = OSC.z)"                                                                       ,
      ""                                                                                        ,
      "model: EM via EM"                                                                        ,
      ""                                                                                        ,
      "              Estimate  l.CI   u.CI"                                                     ,
      "ERR              0.618 0.447  0.743"                                                     ,
      "EDR              0.391 0.070  0.701"                                                     ,
      "Soric FDR        0.082 0.022  0.697"                                                     ,
      "File Drawer R    1.560 0.426 13.237"                                                     ,
      "Expected N         218   121   1210"                                                     ,
      "Missing N          128    31   1120"                                                     ,
      ""                                                                                        ,
      "Model converged in 27 + 783 iterations"                                                  ,
      "Fitted using 73 z-values. 90 supplied, 85 significant (ODR = 0.94, 95% CI [0.87, 0.98]).",
      "Q = -60.61, 95% CI[-72.24, -46.24]"       
    ))
  
  # plot
  skip_on_os("mac")
  vdiffr::expect_doppelganger("z-curve EM", function()plot(fit.EM, main = "OSC (with EM)", annotation = TRUE, CI = TRUE))
  vdiffr::expect_doppelganger("z-curve EM (ggplot)", suppressWarnings(plot(fit.EM, main = "OSC (with EM)", annotation = TRUE, CI = TRUE, plot_type = "ggplot")))

})
test_that("z-curve KD2 can be fitted and reproduces OSC results", {

  # set seed for reproducibility 
  set.seed(666)
  skip_on_os("mac")
  
  # fit the EM method
  fit.KD2 <- suppressWarnings(zcurve(OSC.z, method = "density", bootstrap = 10))

  # Round fitted results to remove small numerical differences across platforms.
  fit.KD2$coefficients <- round(fit.KD2$coefficients, 5)
  fit.KD2$fit$weights <- round(fit.KD2$fit$weights, 5)
  
  # print
  expect_equal(
    capture_output_lines(print(fit.KD2), print = TRUE, width = 150),
    c("Call:"                                                  ,
      "zcurve(z = OSC.z, method = \"density\", bootstrap = 10)",
      ""                                                       ,
      "Estimates:"                                             ,
      "    ERR     EDR "                                         ,
      "0.61327 0.50647 "
    ))
  
  # basic summary
  expect_equal(
    capture_output_lines(summary(fit.KD2), print = TRUE, width = 150),
    c("Call:"                                                                                   ,
      "zcurve(z = OSC.z, method = \"density\", bootstrap = 10)"                                 ,
      ""                                                                                        ,
      "model: KD2 via density"                                                                  ,
      ""                                                                                        ,
      "    Estimate  l.CI  u.CI"                                                                ,
      "ERR    0.613 0.495 0.745"                                                                ,
      "EDR    0.506 0.141 0.714"                                                                ,
      ""                                                                                        ,
      "Model converged in 46 iterations"                                                        ,
      "Fitted using 73 z-values. 90 supplied, 85 significant (ODR = 0.94, 95% CI [0.87, 0.98]).",
      "RMSE = 0.11, 95% CI[0.09, 0.19]"
    ))
  
  # plot
  vdiffr::expect_doppelganger("z-curve KD2", function()plot(fit.KD2))
  vdiffr::expect_doppelganger("z-curve KD2 (ggplot)", plot(fit.KD2, plot_type = "ggplot"))
  
})
test_that("z-curve EM censoring works", {
  
  set.seed(1)
  z    <- runif(100, 2, 7)
  z.lb <- runif(100, 2, 7)
  z.ub <- z.lb + runif(100, 0, 1)
  
  # mixed censored / normal fit
  fit.mixed  <- suppressWarnings(zcurve(z = z, z.lb = z.lb, z.ub = z.ub))
  
  # summary
  expect_equal(
    capture_output_lines(summary(fit.mixed), print = TRUE, width = 150),
    c("Call:"                                                                                      ,
      "zcurve(z = z, z.lb = z.lb, z.ub = z.ub)"                                                    ,
      ""                                                                                           ,
      "model: EM via EM"                                                                           ,
      ""                                                                                           ,
      "    Estimate  l.CI  u.CI"                                                                   ,
      "ERR    0.961 0.910 1.000"                                                                   ,
      "EDR    0.957 0.884 1.000"                                                                   ,
      ""                                                                                           ,
      "Model converged in 70 + 223 iterations"                                                     ,
      "Fitted using 165 z-values. 200 supplied, 200 significant (ODR = 1.00, 95% CI [0.98, 1.00]).",
      "Q = -319.22, 95% CI[-346.77, -288.55]"
    ))
  
  # plot
  vdiffr::expect_doppelganger("z-curve mixed EM", function()plot(fit.mixed, CI = TRUE))
  vdiffr::expect_doppelganger("z-curve mixed EM (ggplot)", plot(fit.mixed, CI = TRUE, plot_type = "ggplot"))
  
  
  # censoring only
  fit.cens  <- suppressWarnings(zcurve(z.lb = z.lb, z.ub = z.ub))
  
  # summary
  expect_equal(
    capture_output_lines(summary(fit.cens), print = TRUE, width = 150),
    c("Call:"                                                                                    ,
      "zcurve(z.lb = z.lb, z.ub = z.ub)"                                                         ,
      ""                                                                                         ,
      "model: EM via EM"                                                                         ,
      ""                                                                                         ,
      "    Estimate  l.CI  u.CI"                                                                 ,
      "ERR    0.968 0.908 1.000"                                                                 ,
      "EDR    0.965 0.879 1.000"                                                                 ,
      ""                                                                                         ,
      "Model converged in 85 + 89 iterations"                                                    ,
      "Fitted using 82 -values. 100 supplied, 100 significant (ODR = 1.00, 95% CI [0.95, 1.00]).",
      "Q = -205.22, 95% CI[-228.49, -179.08]"  
    ))
  
  # plot
  vdiffr::expect_doppelganger("z-curve cens EM", function()plot(fit.cens, CI = TRUE))
  vdiffr::expect_doppelganger("z-curve cens EM (ggplot)", plot(fit.cens, CI = TRUE, plot_type = "ggplot"))
})
test_that("z-curve clustered works", {
  
  set.seed(1)
  z    <- abs(c(rnorm(300, 3, 1), rnorm(200, 1, 1)))
  id   <- sample(500, 500, TRUE)
  
  # rounded
  data <- paste0("z = ",round(z, 2))
  data <- zcurve_data(data, id)
  fitb <- suppressWarnings(zcurve_clustered(data = data, method = "b", bootstrap = 20))
  fitw <- suppressWarnings(zcurve_clustered(data = data, method = "w", bootstrap = 20))
  
  expect_equal(
    capture_output_lines(summary(fitb), print = TRUE, width = 150),
    c("Call:"                                                                                                ,
      "zcurve_clustered(data = data, method = \"b\", bootstrap = 20)"                                        ,
      ""                                                                                                     ,
      "model: EM via EM (bootstrapped)"                                                                      ,
      ""                                                                                                     ,
      "    Estimate  l.CI  u.CI"                                                                             ,
      "ERR    0.760 0.695 0.829"                                                                             ,
      "EDR    0.722 0.632 0.823"                                                                             ,
      ""                                                                                                     ,
      "Model converged in 57 + 450.85 iterations"                                                               ,
      "Fitted using 299 zcurve-data-values. 500 supplied, 300 significant (ODR = 0.60, 95% CI [0.56, 0.64]).",
      "Q = -236.27, 95% CI[-251.36, -224.63]"  
    ))
  
  expect_equal(
    capture_output_lines(summary(fitw), print = TRUE, width = 150),
    c("Call:"                                                                                                ,
      "zcurve_clustered(data = data, method = \"w\", bootstrap = 20)"                                        ,
      ""                                                                                                     ,
      "model: EM via EM (weighted)"                                                                          ,
      ""                                                                                                     ,
      "    Estimate  l.CI  u.CI"                                                                             ,
      "ERR    0.766 0.675 0.861"                                                                             ,
      "EDR    0.730 0.613 0.868"                                                                             ,
      ""                                                                                                     ,
      "Model converged in 73 + 938 iterations"                                                               ,
      "Fitted using 299 zcurve-data-values. 500 supplied, 300 significant (ODR = 0.60, 95% CI [0.56, 0.64]).",
      "Q = -236.38, 95% CI[-255.04, -216.19]"  
    ))
  
  vdiffr::expect_doppelganger("z-curve clustered rounded", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    plot(fitb)
    plot(fitw)
  })

    
  # precise
  data <- paste0("z = ",z)
  data <- zcurve_data(data, id)
  fitb <- suppressWarnings(zcurve_clustered(data = data, method = "b", bootstrap = 20))
  fitw <- suppressWarnings(zcurve_clustered(data = data, method = "w", bootstrap = 20))
  
  expect_equal(
    capture_output_lines(summary(fitb), print = TRUE, width = 150),
    c("Call:"                                                                                                ,
      "zcurve_clustered(data = data, method = \"b\", bootstrap = 20)"                                        ,
      ""                                                                                                     ,
      "model: EM via EM (bootstrapped)"                                                                      ,
      ""                                                                                                     ,
      "    Estimate  l.CI  u.CI"                                                                             ,
      "ERR    0.766 0.699 0.830"                                                                             ,
      "EDR    0.730 0.636 0.824"                                                                             ,
      ""                                                                                                     ,
      "Model converged in 78 + 430.65 iterations"                                                            ,
      "Fitted using 299 zcurve-data-values. 500 supplied, 299 significant (ODR = 0.60, 95% CI [0.55, 0.64]).",
      "Q = -199.11, 95% CI[-209.72, -191.18]"  
    ))
  
  expect_equal(
    capture_output_lines(summary(fitw), print = TRUE, width = 150),
    c("Call:"                                                                                                ,
      "zcurve_clustered(data = data, method = \"w\", bootstrap = 20)"                                        ,
      ""                                                                                                     ,
      "model: EM via EM (weighted)"                                                                          ,
      ""                                                                                                     ,
      "    Estimate  l.CI  u.CI"                                                                             ,
      "ERR    0.766 0.671 0.850"                                                                             ,
      "EDR    0.730 0.610 0.851"                                                                             ,
      ""                                                                                                     ,
      "Model converged in 78 + 676 iterations"                                                               ,
      "Fitted using 299 zcurve-data-values. 500 supplied, 299 significant (ODR = 0.60, 95% CI [0.55, 0.64]).",
      "Q = -199.94, 95% CI[-211.34, -175.66]" 
    ))
  
  vdiffr::expect_doppelganger("z-curve clustered precise", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    plot(fitb)
    plot(fitw)
  })
  
  
  # mixed
  data <- paste0("z = ",c(round(z[1:100], 1), round(z[101:200], 2), z[201:500]))
  data <- zcurve_data(data, id)
  fitb <- suppressWarnings(zcurve_clustered(data = data, method = "b", bootstrap = 20))
  fitw <- suppressWarnings(zcurve_clustered(data = data, method = "w", bootstrap = 20))
  
  expect_equal(
    capture_output_lines(summary(fitb), print = TRUE, width = 150),
    c("Call:"                                                                                                ,
      "zcurve_clustered(data = data, method = \"b\", bootstrap = 20)"                                        ,
      ""                                                                                                     ,
      "model: EM via EM (bootstrapped)"                                                                      ,
      ""                                                                                                     ,
      "    Estimate  l.CI  u.CI"                                                                             ,
      "ERR    0.773 0.712 0.834"                                                                             ,
      "EDR    0.738 0.653 0.830"                                                                             ,
      ""                                                                                                     ,
      "Model converged in 68 + 481.75 iterations"                                                            ,
      "Fitted using 299 zcurve-data-values. 500 supplied, 300 significant (ODR = 0.60, 95% CI [0.56, 0.64]).",
      "Q = -316.94, 95% CI[-328.42, -301.35]"  
    ))
  
  expect_equal(
    capture_output_lines(summary(fitw), print = TRUE, width = 150),
    c("Call:"                                                                                                ,
      "zcurve_clustered(data = data, method = \"w\", bootstrap = 20)"                                        ,
      ""                                                                                                     ,
      "model: EM via EM (weighted)"                                                                          ,
      ""                                                                                                     ,
      "    Estimate  l.CI  u.CI"                                                                             ,
      "ERR    0.768 0.687 0.842"                                                                             ,
      "EDR    0.732 0.618 0.838"                                                                             ,
      ""                                                                                                     ,
      "Model converged in 76 + 274 iterations"                                                               ,
      "Fitted using 299 zcurve-data-values. 500 supplied, 300 significant (ODR = 0.60, 95% CI [0.56, 0.64]).",
      "Q = -319.69, 95% CI[-353.95, -301.27]"  
    ))
  
  vdiffr::expect_doppelganger("z-curve clustered mixed", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    plot(fitb)
    plot(fitw)
  })
  
  vdiffr::expect_doppelganger("z-curve clustered mixed (ggplot-1)", plot(fitb, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("z-curve clustered mixed (ggplot-2)", plot(fitw, plot_type = "ggplot"))
})