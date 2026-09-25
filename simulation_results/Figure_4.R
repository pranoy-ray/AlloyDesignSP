# Figure 4: arrows from SP to the coordinate-scaled, renormalized SSP design.
# Resolve paths for Rscript and source(), including invocation from another directory.
entry_file <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(files)) tail(files, 1)[[1]] else
    sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(entry_file, mustWork = TRUE)), "common.R"), local = TRUE)

require_packages(c(sampling_packages, "Ternary", "magick"))
source(file.path(script_dir, "sampling_strategy.R"), local = TRUE)

p <- 3L
n <- 15L
candidate_set <- unique(simplex_trans(spacefillr::generate_sobol_set(10000, p - 1)))
sp <- get_sp(candidate_set, n, seed = 3)
ssp <- get_ssp(candidate_set, n, seed = 3)
# The common seed preserves row correspondence between SP and SSP.
save_ternary("all_X.png", function() {
  ternary_axes()
  Ternary::TernaryPoints(rbind(sp$X_sub, ssp$X_sub),
                         pch = rep(c(21, 23), each = n),
                         bg = rep(c("red", "green"), each = n), lwd = 0.9, cex = 2)
  for (i in seq_len(n)) {
    Ternary::TernaryArrows(sp$X_sub[i, ], ssp$X_sub[i, ], length = 0.1, col = "darkblue")
  }
  legend("topright", legend = c("Support Points", "Scaled Support Points"),
         cex = 1, bty = "n", pch = c(21, 23), pt.cex = 2, pt.bg = c("red", "green"))
})
saveRDS(list(n = n, seed = 3L, kappa = 1, candidate_set = candidate_set, SP = sp, SSP = ssp),
        file.path(output_dir, "Figure_4_results.rds"))
