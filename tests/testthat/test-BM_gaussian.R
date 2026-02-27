set.seed(12)

test_that("BM_gaussian SBM estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu <- 20 * matrix(runif(Q * Q), Q, Q)
  M <- matrix(rnorm(n * n, sd = 10), n, n) + Z %*% Mu %*% t(Z)

  model <- BM_gaussian("SBM", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_gaussian SBM_sym estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  Mu <- 20 * matrix(runif(Q * Q), Q, Q)
  Mu[lower.tri(Mu)] <- t(Mu)[lower.tri(Mu)]
  M <- matrix(rnorm(n * n, sd = 10), n, n) + Z %*% Mu %*% t(Z)
  M[lower.tri(M)] <- t(M)[lower.tri(M)]

  model <- BM_gaussian("SBM_sym", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_gaussian LBM estimation runs", {
  npc <- c(20, 10)
  Q <- c(1, 2)
  n <- npc * Q
  Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
  Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
  Mu <- 20 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  M <- matrix(rnorm(n[1] * n[2], sd = 10), n[1], n[2]) + Z1 %*% Mu %*% t(Z2)

  model <- BM_gaussian("LBM", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})
