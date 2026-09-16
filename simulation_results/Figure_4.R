# Figure 4: Sensitivity of g over a range of theta
rm(list = ls())
source("sampling_strategy.R")
source("MMSE_calculation.R")

N = 1000
p = 3
n = 6
candidate_set <- MOFAT::qmc_generate(N,p-1) # QMC points in [0,1]^{p-1}
candidate_set <- unique(simplex_trans(candidate_set))

############################################## Support Points
idx = get_sp_idx(candidate_set, n, seed = 1)
sp_X = candidate_set[idx,]

############################################## Scaled Support Points
scaled_SP_g <- function(X, g) {
  a <- lapply(1:ncol(X), function(i) {
    ifelse(X[, i] < 0.5, X[, i] - min(X[, i]) * g,
          ifelse(X[, i] > 0.5, X[, i] + (1 - max(X[, i])) * g, X[, i]))
  })
  b <- Reduce('+', a)
  ans <- lapply(1:ncol(X), function(i) {
    a[[i]] / b
  })
  ans <- Reduce(cbind, ans)
  return(ans)
}

get_ssp_idx_g <- function(X, n, seed, g) {
  set.seed(seed)
  p <- ncol(X)
  X_sp <- support::sp(n, p, dist.samp = X)$sp
  X_scaled_sp <- scaled_SP_g(X_sp, g)
  X_sp_idx <- SPlit::subsample(X, X_scaled_sp) # idx of representative points for X
  return(X_sp_idx)
}

#################### SSP, g = 0
idx = get_ssp_idx_g(candidate_set, n, seed = 1, g = 0)
ssp_X_0 = candidate_set[idx,]

#################### SSP, g = 0.5
idx = get_ssp_idx_g(candidate_set, n, seed = 1, g = 0.5)
ssp_X_05 = candidate_set[idx,]

#################### SSP, g = 1
idx = get_ssp_idx_g(candidate_set, n, seed = 1, g = 1)
ssp_X_1 = candidate_set[idx,]


############################################## MMSE
test_set = MOFAT::qmc_generate(10000,p-1)
test_set <- unique(simplex_trans(test_set))

theta_ls <- seq(1e-1, 2, length.out = 100)
MMSE_ssp_X_0 <- c()
MMSE_ssp_X_05 <- c()
MMSE_ssp_X_1 <- c()
for(i in 1:length(theta_ls)) {
  print(i)
  MMSE_ssp_X_0 <- c(MMSE_ssp_X_0, MMSE(ssp_X_0, test_set, theta_ls[i]))
  MMSE_ssp_X_05 <- c(MMSE_ssp_X_05, MMSE(ssp_X_05, test_set, theta_ls[i]))
  MMSE_ssp_X_1 <- c(MMSE_ssp_X_1, MMSE(ssp_X_1, test_set, theta_ls[i]))
}

############################################## Plot
jpeg("MMSE_sensitivity_of_g.jpeg", width = 7, height = 5,
     units = "in", res = 600, quality = 100)
plot_cols <- c("#0072B2", "#D55E00", "#009E73")
matplot(theta_ls, cbind(MMSE_ssp_X_0, MMSE_ssp_X_05, MMSE_ssp_X_1),
        type = "o", col = plot_cols, lty = 1:3, pch = c(16, 17, 15),
        lwd = 1, cex = 0.4, xlab = expression(theta), ylab = "MMSE")
axis(1, at = theta_ls[1])
legend("topright", legend = c("g = 0", "g = 0.5", "g = 1"),
       col = plot_cols, lty = 1:3, pch = c(16, 17, 15), lwd = 1, bty = "n")
dev.off()