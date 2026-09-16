rm(list = ls())
source("sampling_strategy.R")
source("MMSE_calculation.R")

N = 1000
p = 3
n = 6
candidate_set <- MOFAT::qmc_generate(N,p-1) # QMC points in [0,1]^{p-1}
candidate_set <- unique(simplex_trans(candidate_set))
test_set = MOFAT::qmc_generate(40000,p-1)
test_set <- unique(simplex_trans(test_set))

# Support Points
idx = get_sp_idx(candidate_set, n, seed = 1)
sp_X = candidate_set[idx,]
sp_mse_val <- MSE(sp_X, test_set, 1)

# Support Points
idx = get_ssp_idx(candidate_set, n, seed = 1)
ssp_X = candidate_set[idx,]
ssp_mse_val <- MSE(ssp_X, test_set, 1)


###### Visualization for Support points
png("MSE_ternary_SP.png",width = 7, height = 6, units = "in", res = 600)
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
mse_min <- min(c(sp_mse_val, ssp_mse_val))
mse_max <- max(c(sp_mse_val, ssp_mse_val))
spectrum_bins <- 255
mse_palette <- hcl.colors(spectrum_bins, palette = "viridis")
mse_breaks <- seq(mse_min, mse_max, length.out = length(mse_palette) + 1)
mse_cols <- mse_palette[findInterval(sp_mse_val, mse_breaks, all.inside = TRUE)]
Ternary::TernaryPoints(
  test_set,
  col = mse_cols,
  pch = 16,
  cex = 0.35
)
Ternary::TernaryPoints(
  sp_X,
  pch = 21,
  bg = "red",
  lwd = 0.9,
  cex = 2
)
PlotTools::SpectrumLegend(
  "topright",
  palette = mse_palette,
  legend = round(seq(mse_max, mse_min, length.out = 5), 3),
  title = "MSE",
  bty = "n",
  xpd = NA,
  inset = 0.02,
  cex = 1.5,
  title.cex = 1.7
)
dev.off()



###### Visualization for Scaled Support points
png("MSE_ternary_SSP.png",
    width = 7, height = 6, units = "in", res = 600)
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
mse_min <- min(c(sp_mse_val, ssp_mse_val))
mse_max <- max(c(sp_mse_val, ssp_mse_val))
spectrum_bins <- 255
mse_palette <- hcl.colors(spectrum_bins, palette = "viridis")
mse_breaks <- seq(mse_min, mse_max, length.out = length(mse_palette) + 1)
mse_cols <- mse_palette[findInterval(ssp_mse_val, mse_breaks, all.inside = TRUE)]
Ternary::TernaryPoints(
  test_set,
  col = mse_cols,
  pch = 16,
  cex = 0.35
)
Ternary::TernaryPoints(
  ssp_X,
  pch = 21,
  bg = "red",
  lwd = 0.9,
  cex = 2
)
PlotTools::SpectrumLegend(
  "topright",
  palette = mse_palette,
  legend = round(seq(mse_max, mse_min, length.out = 5), 3),
  title = "MSE",
  bty = "n",
  xpd = NA,
  inset = 0.02,
  cex = 1.5,
  title.cex = 1.7
)
dev.off()