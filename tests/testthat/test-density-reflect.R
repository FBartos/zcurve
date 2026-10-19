# Reference values generated with evmix 2.12:
# dbckden(x, kerncentres, bw = .1, bcmethod = "reflect").
# Keep these fixed so the regression tests do not depend on evmix.
test_that("reflected Gaussian densities reproduce evmix", {
  x <- c(-.1, 0, .05, .1, .3, .5, 1, 1.5, 2)
  centres <- c(0, .05, .05, .4, 1)
  expected <- c(0, 4.41282705662319, 3.97369065961181, 2.90308125806024,
                0.575272775477035, 0.484014735110205, 0.797884572954631,
                2.9734390294686e-06, 0)

  expect_equal(zcurve:::.zcurve_density_reflect(x, centres, .1),
               expected, tolerance = 1e-13)
})

test_that("reflected Gaussian densities preserve evmix's evaluation limits", {
  x <- c(0, .49, .5, 1, 1.5, 1.51)
  expected <- c(0, 0, 1.4867195147343e-05, 3.98942280401433,
                1.4867195147343e-05, 0)

  expect_equal(zcurve:::.zcurve_density_reflect(x, 1, .1),
               expected, tolerance = 1e-13)
})

test_that("reflected Gaussian densities handle non-finite inputs", {
  expect_warning(
    actual <- zcurve:::.zcurve_density_reflect(c(NA_real_, NaN, Inf, -Inf, 0), 0, .1),
    "infinite quantiles set to NA")
  expect_equal(actual, c(NA_real_, NaN, NA_real_, NA_real_, 7.97884560802865),
               tolerance = 1e-13)

  expect_warning(
    actual <- zcurve:::.zcurve_density_reflect(0, c(0, NA, Inf), .1),
    "non-finite kernel centres are dropped")
  expect_equal(actual, 7.97884560802865, tolerance = 1e-13)
})

test_that("reflected Gaussian densities reject invalid centres and bandwidths", {
  expect_error(zcurve:::.zcurve_density_reflect(0, -.1, .1),
               "negative kernel centres not permitted")
  expect_error(zcurve:::.zcurve_density_reflect(0, numeric(), .1),
               "non-empty numeric vector")
  expect_warning(
    expect_error(zcurve:::.zcurve_density_reflect(0, NA_real_, .1),
                 "at least one finite value"),
    "non-finite kernel centres are dropped")

  for(bw in list(NULL, numeric(), 0, -1, NA_real_, Inf, "nrd0", c(.1, .2)))
    expect_error(zcurve:::.zcurve_density_reflect(0, 0, bw), "bw must be")
})

test_that("density_dbc fits reproduce evmix at different lower boundaries", {
  # Frozen results from zcurve(OSC.z, method = "density", bootstrap = FALSE,
  # control = list(density_dbc = TRUE, a = a)) with evmix 2.12.
  references <- list(
    list(a = stats::qnorm(.975),
         weights = c(0, 0.248454993565625, 0.406043716353554,
                     0.345501290080822, 0, 0, 0),
         objective = 0.0837264876023584,
         density_x = c(1.95996398454005, 1.96996398454005, 2.05996398454005,
                       2.45996398454005, 2.95996398454005, 5.99996398454005),
         density_y = c(0.601135715458224, 0.601144373769617, 0.62377420201274,
                       0.53299819913007, 0.562174951252083, 0)),
    list(a = 2.5,
         weights = c(0, 0, 0.455588570606487, 0.544411429393513, 0, 0, 0),
         objective = 0.0764510762577429,
         density_x = c(2.5, 2.51, 2.6, 3, 3.5, 6),
         density_y = c(0.654158214666823, 0.655632956530529, 0.757956078659719,
                       0.994408447880425, 0.489683559469158, 0))
  )

  for(reference in references){
    fit <- suppressWarnings(zcurve(OSC.z, method = "density", bootstrap = FALSE,
                                   control = list(density_dbc = TRUE, a = reference$a)))
    density <- fit$fit$density
    index <- c(1, 2, 11, 51, 101, length(density$x))

    expect_equal(density$x[index], reference$density_x, tolerance = 1e-13)
    expect_equal(density$y[index], reference$density_y, tolerance = 1e-13)
    expect_equal(sum(density$y * diff(density$x)[1]), 1, tolerance = 1e-13)
    expect_equal(fit$fit$weights, reference$weights, tolerance = 1e-5)
    expect_equal(fit$fit$objective, reference$objective, tolerance = 1e-8)
  }
})
