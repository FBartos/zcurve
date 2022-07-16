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
  expect_doppelganger("z-curve EM", function()plot(fit.EM, main = "OSC (with EM)", annotation = TRUE, CI = TRUE))

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
  expect_doppelganger("z-curve KD2", function()plot(fit.KD2))
  
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
      "ERR    0.961 0.915 1.000"                                                                   ,
      "EDR    0.957 0.889 1.000"                                                                   ,
      ""                                                                                           ,
      "Model converged in 70 + 223 iterations"                                                     ,
      "Fitted using 165 z-values. 200 supplied, 200 significant (ODR = 1.00, 95% CI [0.98, 1.00]).",
      "Q = -319.22, 95% CI[-332.77, -303.93]"
    ))
  
  # plot
  expect_doppelganger("z-curve mixed EM", function()plot(fit.mixed, CI = TRUE))
  
  
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
      "ERR    0.968 0.938 0.998"                                                                 ,
      "EDR    0.965 0.915 1.000"                                                                 ,
      ""                                                                                         ,
      "Model converged in 88 + 63 iterations"                                                    ,
      "Fitted using 82 -values. 100 supplied, 100 significant (ODR = 1.00, 95% CI [0.95, 1.00]).",
      "Q = -205.22, 95% CI[-205.22, -205.22]"  
    ))
  
  # plot
  expect_doppelganger("z-curve cens EM", function()plot(fit.cens, CI = TRUE))
})