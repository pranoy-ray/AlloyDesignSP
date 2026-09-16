# MMSE calculation
# Isotropic Gaussian correlation: R(x - y) = exp(-||x - y||^2 / theta^2).
MMSE <- function(D, test_set, theta) {
  stopifnot(length(theta) == 1L, is.finite(theta), theta > 0,
            ncol(D) == ncol(test_set), nrow(D) > 0, nrow(test_set) > 0)
  kernel <- function(A, B) {
    d2 <- outer(rowSums(A^2), rowSums(B^2), "+") - 2 * tcrossprod(A, B)
    exp(-pmax(d2, 0) / theta^2)
  }
  R <- kernel(D, D)
  r <- kernel(D, test_set)
  v <- backsolve(chol(R), r, transpose = TRUE)
  max(0, 1 - colSums(v^2)) # Clamp negative roundoff to zero.
}

# MSE calculation
MSE <- function(D, test_set, theta) {
  stopifnot(length(theta) == 1L, is.finite(theta), theta > 0,
            ncol(D) == ncol(test_set), nrow(D) > 0, nrow(test_set) > 0)
  kernel <- function(A, B) {
    d2 <- outer(rowSums(A^2), rowSums(B^2), "+") - 2 * tcrossprod(A, B)
    exp(-pmax(d2, 0) / theta^2)
  }
  R <- kernel(D, D)
  r <- kernel(D, test_set)
  v <- backsolve(chol(R), r, transpose = TRUE)
  ifelse(1 - colSums(v^2) < 0, 0, 1 - colSums(v^2))
}