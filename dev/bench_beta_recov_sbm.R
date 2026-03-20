library("nnet")
library("blockmodels")
library("aricode")
library("sbm")

n <- 200
p <- 3
K <- 3

alpha_cp <- matrix(c(
  0.9, 0.5, 0.05,
  0.3, 0.2, 0.1,
  0.05, 0.01, 0
), byrow = TRUE, nrow = K)

alpha_as <- matrix(c(
  0.7, 0.1, 0.05,
  0.08, 0.5, 0.02,
  0.05, 0.08, 0.25
), byrow = TRUE, nrow = K)

alpha_dis <- matrix(c(
  0.08, 0.45, 0.9,
  0.8, 0.05, 0.15,
  0.5, 0.7, 0.01
), byrow = TRUE, nrow = K)

alphas_matrices <- list(assortative = alpha_as, core_periphery = alpha_cp, disassortative = alpha_dis)

reorder_beta <- function(fit) {
  max_idx <- which.max(fit$ICL)
  fit$memberships[[max_idx]]$B

  orderLabels <- order(colMeans(fit$memberships[[max_idx]]$alpha) %*%
    fit$model_parameters[[max_idx]]$pi, decreasing = TRUE)
  beta_hat <- matrix(fit$memberships[[max_idx]]$B[, orderLabels], ncol = max_idx)

  beta_tilde <- beta_hat - beta_hat[, ncol(beta_hat)]
  beta_tilde
}

seeds <- seq(1L, 100L)
alphas <- seq_along(alphas_matrices)
conditions <- expand.grid(seed = seeds, alpha = alphas)
library(future.apply)
library(future.callr)
plan("callr")
results_df <- future_lapply(seq_len(nrow(conditions)), function(c_idx) {
  alpha_idx <- conditions[c_idx, "alpha"]
  seed <- conditions[c_idx, "seed"]
  alpha <- alphas_matrices[[alpha_idx]]
  message(c_idx, "/", nrow(conditions), ": ", "seed ", seed, " and ", names(alphas_matrices)[alpha_idx], " structure.")
  set.seed(seed)

  X <- matrix(1, n, p + 1)
  for (j in 2:(p + 1)) {
    X[, j] <- rnorm(n, 0, 1 / 2)
  }

  Beta <- round(matrix(abs(rnorm((p + 1) * K, 0, 3)), (p + 1), K))
  Beta[, K] <- 0
  Sign_Beta <- matrix(0, p + 1, K)
  for (i in 1:(p + 1))
  {
    for (k in 1:K) {
      Sign_Beta[i, k] <- (-1)^(k + i)
    }
  }
  Beta <- Beta * Sign_Beta
  nb_Beta_coeff <- nrow(Beta) * (ncol(Beta) - 1)

  TAU <- exp(X %*% Beta)
  Z <- sapply(1:n, function(i) {
    sample(1:3, size = 1, replace = FALSE, TAU[i, ])
  })

  Y <- matrix(rbinom(n * n, 1, alpha[Z, Z]), n, n)

  # sbm
  res_SBM <- estimateSimpleSBM(Y, model = "bernoulli", estimOptions = list(plot = FALSE, verbosity = 0))
  alpha_sbm <- res_SBM$connectParam$mean
  K_sbm <- ncol(alpha_sbm)
  if (K_sbm != K) {
    return(data.frame(seed = seed, alpha_struct = names(alphas_matrices)[alpha_idx], K_hat = K_sbm, ARI = c(NA, NA, NA), MSE = rep(NA, 3), model = c("sbm_multinom_taus", "sbm_multinom_Z", "blockmodels")))
  }
  ord_sbm <- order(apply(alpha_sbm, 1, sum), decreasing = TRUE)
  tau_sbm <- res_SBM$probMemberships[, ord_sbm]
  Z_sbm <- apply(res_SBM$indMemberships[, ord_sbm], 1, which.max)

  ari_sbm <- ARI(Z, Z_sbm)

  Beta_sbm_multinom_taus <- cbind(0, t(summary(multinom(tau_sbm ~ X[, 2] + X[, 3] + X[, 4], ))$coefficients))
  Beta_sbm_multinom_taus <- Beta_sbm_multinom_taus - Beta_sbm_multinom_taus[, K_sbm]
  Beta_sbm_multinom_Z <- cbind(rep(0, p + 1), t(summary(multinom(Z_sbm ~ X[, 2] + X[, 3] + X[, 4]))$coefficients))
  Beta_sbm_multinom_Z <- Beta_sbm_multinom_Z - Beta_sbm_multinom_Z[, K]

  # blockmodels
  res_bm <- BM_bernoulli(
    membership_type = "SBM", adj = Y,
    plotting = "", nodes_covariates = list(node = X), verbosity = 0
  )
  res_bm$estimate()
  K_bm <- which.max(res_bm$ICL)
  ari_bm <- ARI(Z, res_bm$memberships[[K_bm]]$map()$C)
  Beta_bm <- reorder_beta(res_bm)
  if (K_bm != K) {
    MSE_Betas <- sapply(
      list(Beta_sbm_multinom_taus, Beta_sbm_multinom_Z, matrix(NA, nrow(Beta), ncol(Beta))), function(Beta_test, Beta_ref = Beta, nb = nb_Beta_coeff) {
        (1 / nb) * sum((Beta_ref - Beta_test)^2)
      }
    )
  } else {
    MSE_Betas <- sapply(
      list(Beta_sbm_multinom_taus, Beta_sbm_multinom_Z, Beta_bm), function(Beta_test, Beta_ref = Beta, nb = nb_Beta_coeff) {
        (1 / nb) * sum((Beta_ref - Beta_test)^2)
      }
    )
  }


  return(data.frame(seed = seed, alpha_struct = names(alphas_matrices)[alpha_idx], K_hat = c(K_sbm, K_sbm, K_bm), ARI = c(ari_sbm, ari_sbm, ari_bm), MSE = MSE_Betas, model = c("sbm_multinom_taus", "sbm_multinom_Z", "blockmodels")))
}, future.seed = NULL) |> do.call(what = "rbind")

library(dplyr)
library(tidyr)
library(ggplot2)
plot_data <- results_df %>%
  na.omit() %>%
  filter(ARI == 1)
ggplot(
  plot_data,
  aes(x = as.factor(alpha_struct), y = MSE, fill = as.factor(model))
) +
  geom_violin(position = "dodge")

ggplot(
  plot_data %>% filter(alpha_struct != "disassortative"),
  aes(x = as.factor(alpha_struct), y = MSE, fill = as.factor(model))
) +
  geom_violin(position = "dodge")

# Computing blockmodels vs sbm + multinom on taus

plot_diff_baseline <- plot_data %>% filter(model != "sbm_multinom_Z") %>% pivot_wider(names_from = model, values_from = MSE) %>% mutate(diff_baseline = blockmodels - sbm_multinom_taus) 
ggplot(
  plot_diff_baseline,
  aes(x = as.factor(alpha_struct), y = diff_baseline)
) +
  geom_boxplot()

ggplot(
  plot_diff_baseline %>% filter(alpha_struct != "disassortative"),
  aes(x = as.factor(alpha_struct), y = diff_baseline)
) +
  geom_boxplot()
