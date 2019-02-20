# Load data and get tot_time
data(Cetacea)
options(digits=16)

# Parameters for validation
precision_lh <- 1e-3
precision_aicc <- 5e-2
precision_param <- 5e-2

context("Non constant rate")
result <- fit_coal_cst(Cetacea, tau0=.Machine$double.eps, gamma=-0.5, cst.rate=FALSE, N0=89)

reference_lh <- 25.744
reference_aicc <- -47.341
reference_tau0 <- 0.106
reference_gamma <- 0.112

test_that("Non constant rate",{
  expect_that( abs(result$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(result$aicc - reference_aicc)  , is_less_than(precision_aicc) )
  expect_that( abs(result$tau0 - reference_tau0)  , is_less_than(precision_param) )
  expect_that( abs(result$gamma - reference_gamma)  , is_less_than(precision_param) )
})


context("Variable rate")
result <- fit_coal_var(Cetacea, lamb0=0.01, alpha=-0.001, mu0=0.0, beta=0, N0=89)
reference_lh <- 25.738
reference_aicc <- -42.976
reference_lamb0 <- 0.107
reference_alpha <- 0.0019
reference_mu0 <- 3.3e-07
reference_beta <- -0.0514

test_that("Non constant rate",{
  expect_that( abs(result$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(result$aicc - reference_aicc)  , is_less_than(precision_aicc) )
  expect_that( abs(result$lamb0 - reference_lamb0)  , is_less_than(precision_param) )
  expect_that( abs(result$alpha - reference_alpha)  , is_less_than(precision_param) )
  expect_that( abs(result$mu0 - reference_mu0)  , is_less_than(precision_param) )
  expect_that( abs(result$beta - reference_beta)  , is_less_than(precision_param) )
})
