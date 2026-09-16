# Figure 2: Illustration of different experimental design strategies
rm(list = ls())
source("sampling_strategy.R")

N = 1000
p = 3
n = 6
candidate_set <- MOFAT::qmc_generate(N,p-1) # QMC points in [0,1]^{p-1}
candidate_set <- unique(simplex_trans(candidate_set))

# random sampling
idx = get_rand_idx(candidate_set, n, seed = 12)
rand_X = candidate_set[idx,]

png("rand_X.png", width = 7, height = 6, units = "in", res = 600)
par(mar = c(0.6, 0.6, 0.6, 4.2), xpd = NA)
Ternary::TernaryPlot(
  atip = "x1",
  btip = "x2",
  ctip = "x3",
  lab.cex = 1.7,
  axis.cex = 1.5,
  grid.lines = 5,
  grid.minor.lines = 1,
  grid.lty = "solid",
  grid.minor.lty = "dotted",
  col = NA,
  grid.col = "grey85",
  axis.col = "grey40",
  ticks.col = "grey40",
  axis.labels = seq(0, 1, by = 0.2),
  axis.rotate = FALSE,
  padding = 0.08
)
Ternary::TernaryPoints(
  rand_X,
  pch = 21,
  bg = "red",
  lwd = 0.9,
  cex = 2
)
dev.off()

# QMC sampling
idx = get_sobol_idx(candidate_set, n, seed = 8)
sobol_X = candidate_set[idx,]

png("sobol_X.png", width = 7, height = 6, units = "in", res = 600)
par(mar = c(0.6, 0.6, 0.6, 4.2), xpd = NA)
Ternary::TernaryPlot(
  atip = "x1",
  btip = "x2",
  ctip = "x3",
  lab.cex = 1.7,
  axis.cex = 1.5,
  grid.lines = 5,
  grid.minor.lines = 1,
  grid.lty = "solid",
  grid.minor.lty = "dotted",
  col = NA,
  grid.col = "grey85",
  axis.col = "grey40",
  ticks.col = "grey40",
  axis.labels = seq(0, 1, by = 0.2),
  axis.rotate = FALSE,
  padding = 0.08
)
Ternary::TernaryPoints(
  sobol_X,
  pch = 21,
  bg = "red",
  lwd = 0.9,
  cex = 2
)
dev.off()


# K-medoids clustering
idx = get_k_medoids_idx(candidate_set, n, seed = 8)
k_medoids_X = candidate_set[idx,]

png("k_med_X.png", width = 7, height = 6, units = "in", res = 600)
par(mar = c(0.6, 0.6, 0.6, 4.2), xpd = NA)
Ternary::TernaryPlot(
  atip = "x1",
  btip = "x2",
  ctip = "x3",
  lab.cex = 1.7,
  axis.cex = 1.5,
  grid.lines = 5,
  grid.minor.lines = 1,
  grid.lty = "solid",
  grid.minor.lty = "dotted",
  col = NA,
  grid.col = "grey85",
  axis.col = "grey40",
  ticks.col = "grey40",
  axis.labels = seq(0, 1, by = 0.2),
  axis.rotate = FALSE,
  padding = 0.08
)
Ternary::TernaryPoints(
  k_medoids_X,
  pch = 21,
  bg = "red",
  lwd = 0.9,
  cex = 2
)
dev.off()