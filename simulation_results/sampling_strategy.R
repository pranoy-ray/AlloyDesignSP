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


# Boundary-emphasizing transformation that pushes support points toward the boundary of the simplex.
scaled_SP <- function(X) {
  a <- lapply(1:ncol(X), function(i) {
    ifelse(X[, i] < 0.5, X[, i] - min(X[, i]) / 2,
          ifelse(X[, i] > 0.5, X[, i] + (1 - max(X[, i])) / 2, X[, i]))
  })
  b <- Reduce('+', a)
  ans <- lapply(1:ncol(X), function(i) {
    a[[i]] / b
  })
  ans <- Reduce(cbind, ans)
  return(ans)
}


################################################
# ------------ Random Sampling ------------
get_rand_idx <- function(X, n, seed) {
  set.seed(seed)
  N <- nrow(X)
  rdx <- sample(N, n)
  return(rdx)
}

# ------------ QMC Sampling (Sobol' sequence) ------------
get_sobol_idx <- function(X, n, seed) {
  p <- ncol(X)
  sobol <- spacefillr::generate_sobol_owen_set(n, p-1, seed)
  sobol_simp <- simplex_trans(sobol)
  b <- SPlit::subsample(X, sobol_simp)
  return(b)
}

# ------------ K-medoids clustering ------------
get_k_medoids_idx <- function(X, n, seed) {
  set.seed(seed)
  medoids = cluster::pam(X, n, variant = "faster")
  return(medoids$id.med)
}

# ------------ Support Points ------------
get_sp_idx <- function(X, n, seed) {
  set.seed(seed)
  p <- ncol(X)
  X_sp <- support::sp(n, p, dist.samp = X)$sp
  X_sp_idx <- SPlit::subsample(X, X_sp)
  return(X_sp_idx) # return the row indices of the representative points in X
}

# ------------ Scaled Support Points ------------
get_ssp_idx <- function(X, n, seed) {
  set.seed(seed)
  p <- ncol(X)
  X_sp <- support::sp(n, p, dist.samp = X)$sp
  X_scaled_sp <- scaled_SP(X_sp)
  X_sp_idx <- SPlit::subsample(X, X_scaled_sp) # idx of representative points for X
  return(X_sp_idx)
}