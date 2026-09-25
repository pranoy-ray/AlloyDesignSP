################################################
# Inverse-probability transformation by Wang and Fang (1990)
# The input X should be a uniform design on the (p-1) unit hypercube.
simplex_trans <- function(X) {
  n <- nrow(X)
  p <- ncol(X) + 1
  x_1 <- 1 - X[,1]^(1/(p-1))
  if((p-1) > 1) {
    x_i <- lapply(2:(p-1), function(i) {
      a <- matrix(X[,1:(i-1)], nrow = n)
      a2 <- lapply(1:ncol(a), function(j) {
        a[,j]^(1/(p-j))
      })
      a2 <- matrix(unlist(a2), nrow = n)
      b <- apply(a2, 1, prod)
      b2 <- b*(1 - X[,i]^(1/(p-i)))
      return(b2)
    })
  }else {
    x_i <- NULL
  }
  a <- matrix(X[,1:(p-1)], nrow = n)
  a2 <- lapply(1:ncol(a), function(j) {
    a[,j]^(1/(p-j))
  })
  a2 <- matrix(unlist(a2), nrow = n)
  x_p <- apply(a2, 1, prod)
  
  x_1 <- matrix(x_1, ncol = 1)
  x_p <- matrix(x_p, ncol = 1)
  if(is.null(x_i)) {
    ans <- cbind(x_1, x_p)
  }else {
    ans <- cbind(x_1, Reduce(cbind, x_i), x_p)
  }
  return(ans)
}

################################################
# ------------ Random Sampling ------------
get_rand <- function(X, n, seed) {
  set.seed(seed)
  N <- nrow(X)
  rdx <- sample(N, n)
  ans <- list(idx = rdx, X_sub = X[rdx, , drop = FALSE])
  return(ans)
}

# ------------ QMC Sampling (Sobol' sequence) ------------
get_sobol <- function(X, n, seed) {
  p <- ncol(X)
  sobol <- spacefillr::generate_sobol_set(n, p-1, seed)
  sobol_simp <- simplex_trans(sobol)
  b <- SPlit::subsample(X, sobol_simp)
  ans <- list(idx = b, X_sub = X[b, , drop = FALSE])
  return(ans)
}

# ------------ K-medoids clustering ------------
get_k_medoids <- function(X, n, seed) {
  set.seed(seed)
  medoids <- cluster::pam(X, n, variant = "faster")
  ans <- list(idx = medoids$id.med, X_sub = X[medoids$id.med, , drop = FALSE])
  return(ans)
}

# ------------ Support Points ------------
get_sp <- function(X, n, seed, ini = NA) {
  set.seed(seed)
  p <- ncol(X)
  X_sp <- support::sp(n, p, dist.samp = X, ini = ini)$sp
  X_sp_idx <- SPlit::subsample(X, X_sp)
  list(idx = X_sp_idx, X_sub = X[X_sp_idx, , drop = FALSE])
}


# ------------ Scaled Support Points (paper Eqs. 14-17) ------------
get_ssp <- function(X, n, seed, kappa = 1) {
  if (is.null(kappa)) kappa <- 1
  stopifnot(length(kappa) == 1L, is.finite(kappa), kappa >= 0, kappa <= 1)
  p <- ncol(X)
  sp <- get_sp(X, n, seed)
  if (kappa == 0) return(sp)
  X_sp <- sp$X_sub
  lambda <- vapply(seq_len(p), function(i) {
    a <- X_sp[X_sp[, i] > 1 / p, i]
    b <- X_sp[X_sp[, i] < 1 / p, i]
    bounds <- c((1 - a) / (a - 1 / p), b / (1 / p - b))
    # A coordinate identically at the centroid does not move under scaling.
    if (length(bounds)) min(bounds) else 0
  }, numeric(1)) * kappa
  X_ssp <- sweep(X_sp, 2, 1 + lambda, "*")
  X_ssp <- sweep(X_ssp, 2, lambda / p, "-")
  X_ssp <- pmax(X_ssp, 0) # Keep matrix dimensions while removing boundary roundoff.
  g <- apply(X_ssp,1,sum)
  stopifnot(all(is.finite(g)), all(g > 0))
  X_ssp <- sweep(X_ssp, 1, g, "/")
  X_ssp_idx <- SPlit::subsample(X, X_ssp)
  ans <- list(idx = X_ssp_idx, X_sub = X[X_ssp_idx, , drop = FALSE])
  return(ans)
}
