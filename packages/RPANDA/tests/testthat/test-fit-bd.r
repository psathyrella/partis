# Load data and get tot_time
data(Cetacea)
tot_time<-max(node.age(Cetacea)$ages)
options(digits=17)

# # Parameters for validation
# precision_lh <- 3e-3
# precision_aicc <- 2e-1
# precision_param <- 5e-2

context("B constant")
# B constant
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){0}
lamb_par <- c(0.09)
mu_par <- c()
res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, cst.lamb=TRUE,fix.mu=TRUE, cond="stem", dt=1e-3)

test_that("B constant",{
  # Reference values
  reference_lh <- -279.0280
  reference_aicc <- 560.1031
  reference_lamb <- 0.107
  # Parameters for validation
  precision_lh <- 3e-3
  precision_aicc <- 2e-1
  precision_param <- 5e-2
  expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
  expect_that( abs(res$lamb_par - reference_lamb)  , is_less_than(precision_param) )
})

# BD constant
context("B & D constant")
f.lamb <-function(t,y){y[1]}
f.mu<-function(t,y){y[1]}
lamb_par <- c(0.09)
mu_par <- c(0.3)
res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89,cst.lamb=TRUE,cst.mu=TRUE, cond="stem", dt=1e-3)

test_that("BD constant",{
  # Reference values
  reference_lh <- -279.0280
  reference_aicc <- 562.1989
  reference_lamb <- 0.107
  reference_mu <- -1.0e-7
  # Parameters for validation
  precision_lh <- 3e-3
  precision_aicc <- 2e-1
  precision_param <- 5e-2
  expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
  expect_that( abs(res$lamb_par - reference_lamb)  , is_less_than(precision_param) )
})


# 3) B variable E
context("B exponential")
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
f.mu<-function(t,y){0}
lamb_par<-c(0.05, 0.01)
mu_par<-c()
res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, expo.lamb=TRUE, fix.mu=TRUE, cond="stem", dt=1e-3)

test_that("B variable exponential",{
  # Reference values
  reference_lh <- -278.9887
  reference_aicc <- 562.0485
  reference_lamb <- c(0.103,0.004)
  # Parameters for validation
  precision_lh <- 3e-3
  precision_aicc <- 2e-1
  precision_param <- 5e-2
  expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
  expect_that( abs(res$lamb_par[1] - reference_lamb[1])  , is_less_than(precision_param) )
  expect_that( abs(res$lamb_par[2] - reference_lamb[2])  , is_less_than(precision_param) )
})

# 4) B variable L
context("B linear")
f.lamb <-function(t,y){y[1] + y[2] * t}
f.mu<-function(t,y){0}
lamb_par<-c(0.05, 0.01)
mu_par<-c()
res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, fix.mu=TRUE, cond="stem", dt=1e-3)

test_that("B variable linear",{
  # Reference values
  reference_lh <- -278.9896
  reference_aicc <- 562.0502
  reference_lamb <- c(0.104,0.0004)
  # Parameters for validation
  precision_lh <- 3e-3
  precision_aicc <- 2e-1
  precision_param <- 5e-2
  expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
  expect_that( abs(res$lamb_par[1] - reference_lamb[1])  , is_less_than(precision_param) )
  expect_that( abs(res$lamb_par[2] - reference_lamb[2])  , is_less_than(precision_param) )
})


# 5) B variable E, D constant
context("B exponential, D constant")
f.lamb <-function(t,y){y[1] * exp(y[2] * t)}
mu_par <- c(0.5)
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.05, 0.01)
mu_par<-c(0.1)
res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, expo.lamb=TRUE, cst.mu=TRUE, cond="stem", dt=1e-3)

test_that("B variable exponential, D constant",{
  # Reference values
  reference_lh <- -278.9887
  reference_aicc <- 564.1204
  reference_lamb <- c(0.104,0.004)
  reference_mu <- -1.0e-6
  # Parameters for validation
  precision_lh <- 3e-3
  precision_aicc <- 2e-1
  precision_param <- 5e-2
  expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
  expect_that( abs(res$lamb_par[1] - reference_lamb[1])  , is_less_than(precision_param) )
  expect_that( abs(res$lamb_par[2] - reference_lamb[2])  , is_less_than(precision_param) )
})


# 6) B variable L, D constant
context("B linear, D constant")
f.lamb <-function(t,y){y[1] + y[2] * t}
f.mu<-function(t,y){y[1]}
lamb_par<-c(0.05, 0.01)
mu_par <- c(0.5)
res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, cst.mu= TRUE, cond="stem", dt=1e-3)

test_that("B variable linear, D constant",{
  # Reference values
  reference_lh <- -278.9896
  reference_aicc <- 564.1221
  reference_lamb <- c(0.104,0.0004)
  reference_mu <- 2.5e-7
  # Parameters for validation
  precision_lh <- 3e-3
  precision_aicc <- 2e-1
  precision_param <- 5e-2
  expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
  expect_that( abs(res$lamb_par[1] - reference_lamb[1])  , is_less_than(precision_param) )
  expect_that( abs(res$lamb_par[2] - reference_lamb[2])  , is_less_than(precision_param) )
})


# 7) B constant, D variable E
context("B constant, D exponential")
f.lamb<-function(t,y){y[1]}
f.mu <-function(t,y){y[1] * exp(y[2] * t)}
lamb_par <- c(0.05)
mu_par <-c(0.005, 0.01)
res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, cst.lamb=TRUE, expo.mu=TRUE, cond="stem", dt=1e-3)

test_that("B constant, D exponential",{
  # Reference values
  reference_lh <- -279.0280
  reference_aicc <- 564.1989
  reference_lamb <- 0.107
  reference_mu <- c(-9.637e-09, -0.0111)
  # Parameters for validation
  precision_lh <- 3e-3
  precision_aicc <- 2e-1
  precision_param <- 5e-2
  expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
  expect_that( abs(res$lamb_par - reference_lamb)  , is_less_than(precision_param) )
})


# 8) B constant, D variable L
context("B constant, D linear")
f.lamb<-function(t,y){y[1]}
f.mu <-function(t,y){y[1] + y[2] * t}
lamb_par <- c(0.05)
mu_par <-c(0.005, 0.001)
res <- fit_bd(Cetacea,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=87/89, cst.lamb=TRUE, cond="stem", dt=1e-3)

test_that("B constant, D linear",{
  # Reference values
  reference_lh <- -279.0280
  reference_aicc <- 564.1989
  reference_lamb <- 0.107
  reference_mu <- c(-1.0e-07, 3.2e-08)
  # Parameters for validation
  precision_lh <- 3e-3
  precision_aicc <- 2e-1
  precision_param <- 5e-2
  expect_that( abs(res$LH - reference_lh)  , is_less_than(precision_lh) )
  expect_that( abs(res$aicc - reference_aicc)  , is_less_than(precision_aicc) )
  expect_that( abs(res$lamb_par - reference_lamb)  , is_less_than(precision_param) )
})
