test_that("BM_bernoulli_covariates SBM handles partial NA and is invariant to na_replace_value", {
  set.seed(101)

  npc <- 10
  Q <- 2
  n <- npc * Q
  sigmo <- function(x) 1 / (1 + exp(-x))
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mg <- 8 * matrix(runif(Q * Q), Q, Q) - 4
  Y1 <- matrix(runif(n * n), n, n) - 0.5
  Y2 <- matrix(runif(n * n), n, n) - 0.5
  P <- sigmo(Z %*% Mg %*% t(Z) + 3 * Y1 - 2 * Y2)
  M <- 1 * (matrix(runif(n * n), n, n) < P)

  off_diag <- which(row(M) != col(M))
  set.seed(102)
  M[sample(off_diag, max(1, floor(0.15 * length(off_diag))))] <- NA_real_

  set.seed(103)
  model_default <- BM_bernoulli_covariates("SBM", M, list(Y1, Y2), na_replace_value = 0, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_default$estimate())

  set.seed(103)
  model_custom <- BM_bernoulli_covariates("SBM", M, list(Y1, Y2), na_replace_value = 999, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_custom$estimate())

  expect_equal(model_default$ICL, model_custom$ICL, tolerance = 1e-10)
  k <- which.max(model_default$ICL)
  expect_equal(model_default$model_parameters[[k]]$m, model_custom$model_parameters[[k]]$m, tolerance = 1e-10)
  expect_equal(model_default$model_parameters[[k]]$beta, model_custom$model_parameters[[k]]$beta, tolerance = 1e-10)
})

test_that("BM_bernoulli_covariates_fast SBM handles partial NA and is invariant to na_replace_value", {
  set.seed(111)

  npc <- 10
  Q <- 2
  n <- npc * Q
  sigmo <- function(x) 1 / (1 + exp(-x))
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mg <- 8 * matrix(runif(Q * Q), Q, Q) - 4
  Y1 <- matrix(runif(n * n), n, n) - 0.5
  Y2 <- matrix(runif(n * n), n, n) - 0.5
  P <- sigmo(Z %*% Mg %*% t(Z) + 4 * Y1 - 1.5 * Y2)
  M <- 1 * (matrix(runif(n * n), n, n) < P)

  off_diag <- which(row(M) != col(M))
  set.seed(112)
  M[sample(off_diag, max(1, floor(0.15 * length(off_diag))))] <- NA_real_

  set.seed(113)
  model_default <- BM_bernoulli_covariates_fast("SBM", M, list(Y1, Y2), na_replace_value = 0, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_default$estimate())

  set.seed(113)
  model_custom <- BM_bernoulli_covariates_fast("SBM", M, list(Y1, Y2), na_replace_value = -321, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_custom$estimate())

  expect_equal(model_default$ICL, model_custom$ICL, tolerance = 1e-10)
  k <- which.max(model_default$ICL)
  expect_equal(model_default$model_parameters[[k]]$m, model_custom$model_parameters[[k]]$m, tolerance = 1e-10)
  expect_equal(model_default$model_parameters[[k]]$beta, model_custom$model_parameters[[k]]$beta, tolerance = 1e-10)
})

test_that("BM_poisson_covariates SBM handles partial NA and is invariant to na_replace_value", {
  set.seed(121)

  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  L <- 50 * matrix(runif(Q * Q), Q, Q)
  Y1 <- matrix(runif(n * n), n, n)
  Y2 <- matrix(runif(n * n), n, n)
  E <- Z %*% L %*% t(Z) * exp(2.2 * Y1 - 0.8 * Y2)
  M <- matrix(rpois(length(as.vector(E)), as.vector(E)), n, n)

  off_diag <- which(row(M) != col(M))
  set.seed(122)
  M[sample(off_diag, max(1, floor(0.15 * length(off_diag))))] <- NA_real_

  set.seed(123)
  model_default <- BM_poisson_covariates("SBM", M, list(Y1, Y2), na_replace_value = 0, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_default$estimate())

  set.seed(123)
  model_custom <- BM_poisson_covariates("SBM", M, list(Y1, Y2), na_replace_value = 1234, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_custom$estimate())

  expect_equal(model_default$ICL, model_custom$ICL, tolerance = 1e-10)
  k <- which.max(model_default$ICL)
  expect_equal(model_default$model_parameters[[k]]$lambda, model_custom$model_parameters[[k]]$lambda, tolerance = 1e-10)
  expect_equal(model_default$model_parameters[[k]]$beta, model_custom$model_parameters[[k]]$beta, tolerance = 1e-10)
})

test_that("BM_gaussian_covariates SBM handles partial NA and is invariant to na_replace_value", {
  set.seed(131)

  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu <- 15 * matrix(runif(Q * Q), Q, Q)
  Y1 <- matrix(runif(n * n), n, n)
  M <- matrix(rnorm(n * n, sd = 5), n, n) + Z %*% Mu %*% t(Z) + 2.2 * Y1

  off_diag <- which(row(M) != col(M))
  set.seed(132)
  M[sample(off_diag, max(1, floor(0.05 * length(off_diag))))] <- NA_real_

  set.seed(133)
  model_default <- BM_gaussian_covariates("SBM", M, list(Y1), na_replace_value = 0, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_default$estimate())

  set.seed(133)
  model_custom <- BM_gaussian_covariates("SBM", M, list(Y1), na_replace_value = -123, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_custom$estimate())

  expect_equal(model_default$ICL, model_custom$ICL, tolerance = 1e-10)
  k <- which.max(model_default$ICL)
  expect_equal(model_default$model_parameters[[k]]$mu, model_custom$model_parameters[[k]]$mu, tolerance = 1e-10)
  expect_equal(model_default$model_parameters[[k]]$beta, model_custom$model_parameters[[k]]$beta, tolerance = 1e-10)
  expect_equal(model_default$model_parameters[[k]]$sigma2, model_custom$model_parameters[[k]]$sigma2, tolerance = 1e-10)
})

test_that("BM_bernoulli_multiplex SBM handles partial NA and is invariant to na_replace_value", {
  set.seed(141)

  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)

  P00 <- matrix(runif(Q * Q), Q, Q)
  P10 <- matrix(runif(Q * Q), Q, Q)
  P01 <- matrix(runif(Q * Q), Q, Q)
  P11 <- matrix(runif(Q * Q), Q, Q)
  S <- P00 + P10 + P01 + P11
  P00 <- P00 / S
  P01 <- P01 / S
  P10 <- P10 / S
  P11 <- P11 / S

  U <- matrix(runif(n * n), n, n)
  M1 <- 1 * (U > Z %*% (P00 + P01) %*% t(Z))
  M2 <- 1 * ((U > Z %*% P00 %*% t(Z)) & (U < Z %*% (P00 + P01 + P11) %*% t(Z)))

  off_diag <- which(row(M1) != col(M1))
  set.seed(142)
  idx <- sample(off_diag, max(1, floor(0.15 * length(off_diag))))
  M1[idx] <- NA_real_
  M2[idx] <- NA_real_

  set.seed(143)
  model_default <- BM_bernoulli_multiplex("SBM", list(M1, M2), na_replace_value = 0, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_default$estimate())

  set.seed(143)
  model_custom <- BM_bernoulli_multiplex("SBM", list(M1, M2), na_replace_value = 17, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_custom$estimate())

  expect_equal(model_default$ICL, model_custom$ICL, tolerance = 1e-10)
  k <- which.max(model_default$ICL)
  expect_equal(model_default$model_parameters[[k]]$pi, model_custom$model_parameters[[k]]$pi, tolerance = 1e-10)
})

test_that("BM_gaussian_multivariate SBM handles partial NA and is invariant to na_replace_value", {
  set.seed(151)

  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu1 <- 4 * matrix(runif(Q * Q), Q, Q)
  Mu2 <- 4 * matrix(runif(Q * Q), Q, Q)
  N1 <- matrix(rnorm(n * n, sd = 1), n, n)
  N2 <- matrix(rnorm(n * n, sd = 1), n, n)
  M1 <- Z %*% Mu1 %*% t(Z) + N1
  M2 <- Z %*% Mu2 %*% t(Z) + 8 * N1 + N2

  off_diag <- which(row(M1) != col(M1))
  set.seed(152)
  idx <- sample(off_diag, max(1, floor(0.15 * length(off_diag))))
  M1[idx] <- NA_real_
  M2[idx] <- NA_real_

  set.seed(153)
  model_default <- BM_gaussian_multivariate("SBM", list(M1, M2), na_replace_value = 0, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_default$estimate())

  set.seed(153)
  model_custom <- BM_gaussian_multivariate("SBM", list(M1, M2), na_replace_value = -456, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_custom$estimate())

  expect_equal(model_default$ICL, model_custom$ICL, tolerance = 1e-10)
  k <- which.max(model_default$ICL)
  expect_equal(model_default$model_parameters[[k]]$mu, model_custom$model_parameters[[k]]$mu, tolerance = 1e-10)
  expect_equal(model_default$model_parameters[[k]]$Sigma, model_custom$model_parameters[[k]]$Sigma, tolerance = 1e-10)
})

test_that("BM_gaussian_multivariate_independent SBM handles partial NA and is invariant to na_replace_value", {
  set.seed(161)

  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu1 <- 6 * matrix(runif(Q * Q), Q, Q)
  Mu2 <- 6 * matrix(runif(Q * Q), Q, Q)
  M1 <- matrix(rnorm(n * n, sd = 4), n, n) + Z %*% Mu1 %*% t(Z)
  M2 <- matrix(rnorm(n * n, sd = 7), n, n) + Z %*% Mu2 %*% t(Z)

  off_diag <- which(row(M1) != col(M1))
  set.seed(162)
  idx <- sample(off_diag, max(1, floor(0.15 * length(off_diag))))
  M1[idx] <- NA_real_
  M2[idx] <- NA_real_

  set.seed(163)
  model_default <- BM_gaussian_multivariate_independent("SBM", list(M1, M2), na_replace_value = 0, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_default$estimate())

  set.seed(163)
  model_custom <- BM_gaussian_multivariate_independent("SBM", list(M1, M2), na_replace_value = 73, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_custom$estimate())

  expect_equal(model_default$ICL, model_custom$ICL, tolerance = 1e-10)
  k <- which.max(model_default$ICL)
  expect_equal(model_default$model_parameters[[k]]$mu, model_custom$model_parameters[[k]]$mu, tolerance = 1e-10)
  expect_equal(model_default$model_parameters[[k]]$sigma2, model_custom$model_parameters[[k]]$sigma2, tolerance = 1e-10)
})

test_that("BM_gaussian_multivariate_independent_homoscedastic SBM handles partial NA and is invariant to na_replace_value", {
  set.seed(171)

  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu1 <- 5 * matrix(runif(Q * Q), Q, Q)
  Mu2 <- 5 * matrix(runif(Q * Q), Q, Q)
  M1 <- matrix(rnorm(n * n, sd = 4), n, n) + Z %*% Mu1 %*% t(Z)
  M2 <- matrix(rnorm(n * n, sd = 4), n, n) + Z %*% Mu2 %*% t(Z)

  off_diag <- which(row(M1) != col(M1))
  set.seed(172)
  idx <- sample(off_diag, max(1, floor(0.15 * length(off_diag))))
  M1[idx] <- NA_real_
  M2[idx] <- NA_real_

  set.seed(173)
  model_default <- BM_gaussian_multivariate_independent_homoscedastic("SBM", list(M1, M2), na_replace_value = 0, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_default$estimate())

  set.seed(173)
  model_custom <- BM_gaussian_multivariate_independent_homoscedastic("SBM", list(M1, M2), na_replace_value = -1e6, plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_no_error(model_custom$estimate())

  expect_equal(model_default$ICL, model_custom$ICL, tolerance = 1e-10)
  k <- which.max(model_default$ICL)
  expect_equal(model_default$model_parameters[[k]]$mu, model_custom$model_parameters[[k]]$mu, tolerance = 1e-10)
  expect_equal(model_default$model_parameters[[k]]$sigma2, model_custom$model_parameters[[k]]$sigma2, tolerance = 1e-10)
})

test_that("Gaussian multivariate SBM models fail gracefully when all values are NA", {
  n <- 20
  M1 <- matrix(NA_real_, n, n)
  M2 <- matrix(NA_real_, n, n)

  model_mv <- BM_gaussian_multivariate("SBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_error(model_mv$estimate(), "No valid non-missing off-diagonal adjacency values in gaussian_multivariate SBM network\\.")

  model_mvi <- BM_gaussian_multivariate_independent("SBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_error(model_mvi$estimate(), "No valid non-missing off-diagonal adjacency values in gaussian_multivariate_independent SBM network\\.")

  model_mvih <- BM_gaussian_multivariate_independent_homoscedastic("SBM", list(M1, M2), plotting = "", explore_min = 2, explore_max = 2, ncores = 1, verbosity = 0)
  expect_error(model_mvih$estimate(), "No valid non-missing off-diagonal adjacency values in gaussian_multivariate_independent_homoscedastic SBM network\\.")
})