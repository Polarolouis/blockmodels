set.seed(12)

test_that("BM_bernoulli SBM estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  P <- matrix(runif(Q * Q), Q, Q)
  M <- 1 * (matrix(runif(n * n), n, n) < Z %*% P %*% t(Z))

  model <- BM_bernoulli("SBM", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_bernoulli SBM_sym estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  P <- matrix(runif(Q * Q), Q, Q)
  P[lower.tri(P)] <- t(P)[lower.tri(P)]
  M <- 1 * (matrix(runif(n * n), n, n) < Z %*% P %*% t(Z))
  M[lower.tri(M)] <- t(M)[lower.tri(M)]

  model <- BM_bernoulli("SBM_sym", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_bernoulli LBM estimation runs", {
  npc <- c(20, 10)
  Q <- c(1, 2)
  n <- npc * Q
  Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
  Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
  P <- matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  M <- 1 * (matrix(runif(n[1] * n[2]), n[1], n[2]) < Z1 %*% P %*% t(Z2))

  model <- BM_bernoulli("LBM", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})
