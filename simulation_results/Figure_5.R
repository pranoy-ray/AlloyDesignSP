# Figure 5: Visualization of SP and SSP
rm(list = ls())
source("sampling_strategy.R")

N = 1000
p = 3
n = 6
candidate_set <- MOFAT::qmc_generate(N,p-1) # QMC points in [0,1]^{p-1}
candidate_set <- unique(simplex_trans(candidate_set))

# Support Points
idx = get_sp_idx(candidate_set, n, seed = 1)
sp_X = candidate_set[idx,]


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
  sp_X,
  pch = 21,
  bg = "red",
  lwd = 0.9,
  cex = 2
)

# Scaled Support Points
idx = get_ssp_idx(candidate_set, n, seed = 1)
ssp_X = candidate_set[idx,]

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
  ssp_X,
  pch = 21,
  bg = "red",
  lwd = 0.9,
  cex = 2
)

##################################### 
all_X = rbind(sp_X, ssp_X)
png("all_X.png", width = 7, height = 6, units = "in", res = 600)
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
point_cols <- c(rep("red", 6), rep("green", 6))
point_pchs <- c(rep(21, 6), rep(23, 6))
Ternary::TernaryPoints(
  all_X,
  pch = point_pchs,
  bg = point_cols,
  lwd = 0.9,
  cex = 2
)
# Add an arrow
Ternary::TernaryArrows(all_X[1,], all_X[7,], length = 0.1, col = "darkblue")
Ternary::TernaryArrows(all_X[2,], all_X[8,], length = 0.1, col = "darkblue")
Ternary::TernaryArrows(all_X[3,], all_X[9,], length = 0.1, col = "darkblue")
Ternary::TernaryArrows(all_X[4,], all_X[10,], length = 0.1, col = "darkblue")
Ternary::TernaryArrows(all_X[5,], all_X[11,], length = 0.1, col = "darkblue")
Ternary::TernaryArrows(all_X[6,], all_X[12,], length = 0.1, col = "darkblue")

legend("topright", 
       legend = c("Support Points", "Scaled Support Points"),
       cex = 1, bty = "n", pch = c(21,23), pt.cex = 2,
       pt.bg = c("red", "green"),
       )
dev.off()