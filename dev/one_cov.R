library(nnet)

npc <- 100
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

indata <- data.frame(Z = as.factor(Z_factor), X)
indata$Z <- relevel(indata$Z, ref = paste(ncol(B)))

fit_multinom <- multinom(Z ~ 0+X, data=indata)
summary(fit_multinom)

P <- matrix(runif(Q * Q), Q, Q)
M <- 1 * (matrix(runif(n * n), n, n) < Z %*% P %*% t(Z)) ## adjacency matrix

devtools::load_all()
fit <- BM_bernoulli(membership_type = "SBM", adj = M, plotting = character(0), ncores = 1L, nodes_covariates = list(node = matrix(X, ncol = 1)), verbosity = 6)
fit$estimate()
fit$memberships[[which.max(fit$ICL)]]$B
