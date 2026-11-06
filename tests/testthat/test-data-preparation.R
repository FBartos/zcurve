context("data preparation")

test_that("decimal places are read from the reported value", {
  expect_equal(zcurve:::.n_decimals(c("0.03", "2", "2.10", ".05", "1.5e-5", "1e-12", "1e3")),
               c(2, 0, 2, 2, 6, 12, 0))
})

test_that("rounded censored test statistics use their own lower bound", {

  # 'F(2,23) > 10' is rounded to 0 decimals -> F > 9.5
  alone <- zcurve_data("F(2,23) > 10")
  after <- zcurve_data(c("z = 2.1", "F(2,23) > 10"))

  expect_equal(alone$censored$p.ub, stats::pf(9.5, 2, 23, lower.tail = FALSE))
  expect_equal(alone$censored$p.rep, stats::pf(10, 2, 23, lower.tail = FALSE))
  expect_equal(after$censored$p.ub[2], alone$censored$p.ub)

  expect_equal(zcurve_data("t(30) > 3")$censored$p.ub, stats::pt(2.5, 30, lower.tail = FALSE) * 2)
  expect_equal(zcurve_data("z > 4")$censored$p.ub,     stats::pnorm(3.5, lower.tail = FALSE) * 2)
})

test_that("values in scientific notation are not treated as heavily rounded", {

  precise <- zcurve_data("p = 1e-5")
  expect_equal(nrow(precise$precise), 1)
  expect_equal(precise$precise$p, 1e-5)

  censored <- zcurve_data("p < 1e-12")
  expect_equal(censored$censored$p.ub, 1e-12)
})

test_that("bootstrap samples are never split over more cores than samples", {

  old_cores <- zcurve.get_option("max_cores")
  on.exit(zcurve.options(max_cores = old_cores))
  zcurve.options(max_cores = 4)

  expect_equal(zcurve:::.bootstrap_core_load(2),  c(1, 1))
  expect_equal(zcurve:::.bootstrap_core_load(10), c(3, 3, 2, 2))
  expect_equal(sum(zcurve:::.bootstrap_core_load(1001)), 1001)
})

test_that("significant results are counted at sig_level, not at the fitting range", {

  set.seed(1)
  z   <- abs(rnorm(500, 2.5))
  fit <- suppressWarnings(zcurve(z, bootstrap = FALSE, control = list(a = 2.5)))

  expect_equal(fit$N_sig, sum(z >= stats::qnorm(.975)))
})

test_that("zcurve_clustered keeps censored results above the fitting range", {

  skip_on_cran()

  set.seed(1)
  z   <- abs(rnorm(200, 2.5))
  dat <- zcurve_data(c(sprintf("z = %.3f", z), rep("z > 7.00", 50)), id = 1:250)

  fit_z <- suppressWarnings(zcurve(data = dat, bootstrap = FALSE))
  fit_c <- suppressWarnings(zcurve_clustered(dat, method = "w", bootstrap = FALSE))

  expect_equal(length(fit_c$data), 250)
  expect_equal(length(fit_c$data_id), 250)
  expect_equal(fit_c$fit$prop_high, fit_z$fit$prop_high)
})
