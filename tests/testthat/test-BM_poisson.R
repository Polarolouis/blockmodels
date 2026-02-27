set.seed(12)

test_that("BM_poisson SBM estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  L <- 70 * matrix(runif(Q * Q), Q, Q)
  M_in_expectation <- Z %*% L %*% t(Z)
  M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n, n)

  model <- BM_poisson("SBM", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_poisson SBM_sym estimation runs", {
  npc <- 10
  Q <- 2
  n <- npc * Q
  Z <- diag(Q) %x% matrix(1, npc, 1)
  L <- 70 * matrix(runif(Q * Q), Q, Q)
  L[lower.tri(L)] <- t(L)[lower.tri(L)]
  M_in_expectation <- Z %*% L %*% t(Z)
  M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n, n)
  M[lower.tri(M)] <- t(M)[lower.tri(M)]

  model <- BM_poisson("SBM_sym", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})

test_that("BM_poisson LBM estimation runs", {
  npc <- c(20, 10)
  Q <- c(1, 2)
  n <- npc * Q
  Z1 <- diag(Q[1]) %x% matrix(1, npc[1], 1)
  Z2 <- diag(Q[2]) %x% matrix(1, npc[2], 1)
  L <- 70 * matrix(runif(Q[1] * Q[2]), Q[1], Q[2])
  M_in_expectation <- Z1 %*% L %*% t(Z2)
  M <- matrix(rpois(length(as.vector(M_in_expectation)), as.vector(M_in_expectation)), n[1], n[2])

  model <- BM_poisson("LBM", M, plotting = "", explore_min = 2, explore_max = 2, ncores = 2, verbosity = 0)
  expect_model_estimation(model)
})
