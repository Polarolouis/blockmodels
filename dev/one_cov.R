library(nnet)

npc <- 150
Q <- 3
n <- npc * Q

B <- matrix(c(-1, 1, 0), byrow = TRUE, nrow = 1)

set.seed(123);X <- matrix(rnorm(n), ncol = 1)

softmax_rows <- function(x){
    x_shift <- x - apply(x, 1, max)
    ex <- exp(x_shift)
    ex / rowSums(ex)
}

pZ <- softmax_rows(X %*%B)
apply(pZ, 1, which.max)

Z_factor <- sapply(seq_len(nrow(pZ)), function(idx) {
    sample(x = seq_len(ncol(pZ)), size = 1, prob = pZ[idx, ])
})

Z <- t(sapply(seq_along(Z_factor), function(idx) {
    vec <- rep(0, ncol(pZ))
    vec[Z_factor[idx]] <- 1
    vec
}))

# Exp with Sophie's discovery
taus <- pZ  # + rnorm(npc * ncol(pZ), sd = 3) taus <- softmax_rows(taus)

unscaled_taus <- taus/taus[,ncol(taus)]

(t(X)%*%X)^(-1)%*%t(X) %*% log(unscaled_taus)


indata <- data.frame(Z = as.factor(Z_factor), X)
indata$Z <- relevel(indata$Z, ref = paste(ncol(B)))

fit_multinom <- multinom(Z ~ 0+X, data=indata)
summary(fit_multinom)

P <- matrix(c(0.9, 0.5, 0.05,
              0.3, 0.2, 0.1,
              0.05, 0.01, 0), byrow = TRUE, nrow = Q)
M <- 1 * (matrix(runif(n * n), n, n) < Z %*% P %*% t(Z)) ## adjacency matrix

devtools::load_all()
library(microbenchmark)
mb_explicit <- microbenchmark("explicit" = {
fit <- BM_bernoulli(membership_type = "SBM", adj = M, plotting = '', ncores = 1L, nodes_covariates = list(node = matrix(X, ncol = 1)), verbosity = 0)
fit$estimate()
}, times = 10L)

devtools::load_all()
mb_optim <- microbenchmark("optim" = {
fit_optim <- BM_bernoulli(membership_type = "SBM", adj = M, plotting = '', ncores = 1L, nodes_covariates = list(node = matrix(X, ncol = 1)), verbosity = 1)
fit_optim$estimate()
}, times = 10L)

reorder_beta <- function(fit) {
max_idx <- which.max(fit$ICL)
fit$memberships[[max_idx]]$B

orderLabels <- order(colMeans(fit$memberships[[max_idx]]$alpha) %*%
fit$model_parameters[[max_idx]]$pi, decreasing = TRUE)
beta_hat <- matrix(fit$memberships[[max_idx]]$B[,orderLabels],ncol = max_idx)

beta_tilde <- beta_hat - beta_hat[,ncol(beta_hat)]
beta_tilde
}

reorder_beta(fit_optim)
