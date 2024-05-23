context("z-curve")
skip_on_cran()

test_that("z-curve EM can be fitted and reproduces OSC results", {

  # set seed for reproducibility 
  set.seed(666)
  
  # fit the EM method
  fit.EM  <- zcurve(OSC.z)

  # print
  expect_equal(
    capture_output_lines(print(fit.EM), print = TRUE, width = 150),
    c("Call:"               ,
      "zcurve(z = OSC.z)"   ,
      ""                    , 
      "Estimates:"          , 
      "      ERR       EDR ",
      "0.6145784 0.3884366 "    
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
      "ERR    0.615 0.443 0.740"                                                                ,
      "EDR    0.388 0.070 0.699"                                                                ,
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
      "ERR              0.615 0.443  0.740"                                                     ,
      "EDR              0.388 0.070  0.699"                                                     ,
      "Soric FDR        0.083 0.023  0.705"                                                     ,
      "File Drawer R    1.574 0.430 13.387"                                                     ,
      "Expected N         219   122   1223"                                                     ,
      "Missing N          129    32   1133"                                                     ,
      ""                                                                                        ,
      "Model converged in 27 + 783 iterations"                                                  ,
      "Fitted using 73 z-values. 90 supplied, 85 significant (ODR = 0.94, 95% CI [0.87, 0.98]).",
      "Q = -60.61, 95% CI[-72.24, -46.24]"       
    ))
  
  # plot
  vdiffr::expect_doppelganger("z-curve EM", function()plot(fit.EM, main = "OSC (with EM)", annotation = TRUE, CI = TRUE))
  vdiffr::expect_doppelganger("z-curve EM (ggplot)", suppressWarnings(plot(fit.EM, main = "OSC (with EM)", annotation = TRUE, CI = TRUE, plot_type = "ggplot")))

})
test_that("z-curve KD2 can be fitted and reproduces OSC results", {
  
  # set seed for reproducibility 
  set.seed(666)
  
  # fit the EM method
  fit.KD2 <- zcurve(OSC.z, method = "density", bootstrap = 10)
  
  # print
  expect_equal(
    capture_output_lines(print(fit.KD2), print = TRUE, width = 150),
    c("Call:"                                                  ,
      "zcurve(z = OSC.z, method = \"density\", bootstrap = 10)",
      ""                                                       ,
      "Estimates:"                                             ,
      "      ERR       EDR "                                   ,
      "0.6133950 0.5064392 "    
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
      "ERR    0.613 0.496 0.745"                                                                ,
      "EDR    0.506 0.141 0.714"                                                                ,
      ""                                                                                        ,
      "Model converged in 47 iterations"                                                        ,
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
  fit.mixed  <- zcurve(z = z, z.lb = z.lb, z.ub = z.ub)
  
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
  fit.cens  <- zcurve(z.lb = z.lb, z.ub = z.ub)
  
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
      "ERR    0.806 0.714 0.877"                                                                             ,
      "EDR    0.783 0.655 0.895"                                                                             ,
      ""                                                                                                     ,
      "Model converged in 73 + 178 iterations"                                                               ,
      "Fitted using 299 zcurve-data-values. 500 supplied, 300 significant (ODR = 0.60, 95% CI [0.56, 0.64]).",
      "Q = -238.11, 95% CI[-269.95, -219.12]"  
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
      "ERR    0.767 0.691 0.831"                                                                             ,
      "EDR    0.730 0.627 0.825"                                                                             ,
      ""                                                                                                     ,
      "Model converged in 66 + 629.95 iterations"                                                            ,
      "Fitted using 299 zcurve-data-values. 500 supplied, 299 significant (ODR = 0.60, 95% CI [0.55, 0.64]).",
      "Q = -199.02, 95% CI[-209.72, -191.18]"  
    ))
  
  expect_equal(
    capture_output_lines(summary(fitw), print = TRUE, width = 150),
    c("Call:"                                                                                                ,
      "zcurve_clustered(data = data, method = \"w\", bootstrap = 20)"                                        ,
      ""                                                                                                     ,
      "model: EM via EM (weighted)"                                                                          ,
      ""                                                                                                     ,
      "    Estimate  l.CI  u.CI"                                                                             ,
      "ERR    0.809 0.678 0.876"                                                                             ,
      "EDR    0.787 0.617 0.893"                                                                             ,
      ""                                                                                                     ,
      "Model converged in 73 + 177 iterations"                                                               ,
      "Fitted using 299 zcurve-data-values. 500 supplied, 299 significant (ODR = 0.60, 95% CI [0.55, 0.64]).",
      "Q = -201.69, 95% CI[-219.18, -185.82]" 
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
      "ERR    0.769 0.703 0.835"                                                                             ,
      "EDR    0.733 0.639 0.831"                                                                             ,
      ""                                                                                                     ,
      "Model converged in 72 + 430.1 iterations"                                                             ,
      "Fitted using 299 zcurve-data-values. 500 supplied, 300 significant (ODR = 0.60, 95% CI [0.56, 0.64]).",
      "Q = -319.14, 95% CI[-340.91, -306.78]"  
    ))
  
  expect_equal(
    capture_output_lines(summary(fitw), print = TRUE, width = 150),
    c("Call:"                                                                                                ,
      "zcurve_clustered(data = data, method = \"w\", bootstrap = 20)"                                        ,
      ""                                                                                                     ,
      "model: EM via EM (weighted)"                                                                          ,
      ""                                                                                                     ,
      "    Estimate  l.CI  u.CI"                                                                             ,
      "ERR    0.810 0.709 0.881"                                                                             ,
      "EDR    0.788 0.651 0.901"                                                                             ,
      ""                                                                                                     ,
      "Model converged in 59 + 188 iterations"                                                               ,
      "Fitted using 299 zcurve-data-values. 500 supplied, 300 significant (ODR = 0.60, 95% CI [0.56, 0.64]).",
      "Q = -321.39, 95% CI[-361.78, -298.20]"  
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