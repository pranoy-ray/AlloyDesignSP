# Figure 3: MSE surfaces for SP and SSP, with a shared color scale.
# Resolve paths for Rscript and source(), including invocation from another directory.
entry_file <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(files)) tail(files, 1)[[1]] else
    sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(entry_file, mustWork = TRUE)), "common.R"), local = TRUE)

require_packages(c(sampling_packages, ternary_packages))
source(file.path(script_dir, "sampling_strategy.R"), local = TRUE)
source(file.path(script_dir, "MMSE_calculation.R"), local = TRUE)

p <- 3L
n <- 15L # Matches the plotted panels; the paper caption instead says n = 12.
candidate_set <- unique(simplex_trans(spacefillr::generate_sobol_set(10000, p - 1)))
# Evaluate Sobol points on the simplex and explicitly include its vertices.
test_set <- simplex_trans(spacefillr::generate_sobol_set(40000, p - 1))
test_set <- rbind(test_set, diag(p))
theta <- 0.45
designs <- list(SP = get_sp(candidate_set, n, seed = 3),
                SSP = get_ssp(candidate_set, n, seed = 3, kappa = 1))
mse <- lapply(designs, function(design) MSE(design$X_sub, test_set, theta))
mse_range <- range(unlist(mse))
mse_palette <- hcl.colors(255, palette = "YlOrRd", rev = TRUE)
mse_breaks <- seq(mse_range[1], mse_range[2], length.out = length(mse_palette) + 1)
filenames <- c(SP = "MSE_ternary_SP.png", SSP = "MSE_ternary_SSP_kappa_1.png")
for (method in names(designs)) {
  save_ternary(filenames[[method]], function() {
    ternary_axes()
    cols <- mse_palette[findInterval(mse[[method]], mse_breaks, all.inside = TRUE)]
    Ternary::TernaryPoints(test_set, col = cols, pch = 16, cex = 0.35)
    Ternary::TernaryPoints(designs[[method]]$X_sub, pch = 21, bg = "green", lwd = 0.9, cex = 2)
    PlotTools::SpectrumLegend(
      "topright", palette = mse_palette,
      legend = round(seq(mse_range[2], mse_range[1], length.out = 5), 3),
      title = "MSE", bty = "n", xpd = NA, inset = 0.02, cex = 1.5, title.cex = 1.7
    )
  })
}
saveRDS(list(n = n, seed = 3L, theta = theta, kappa = 1, candidate_set = candidate_set,
             test_set = test_set, designs = designs, mse = mse, color_limits = mse_range),
        file.path(output_dir, "Figure_3_results.rds"))
