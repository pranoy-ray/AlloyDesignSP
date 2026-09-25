# Figure 2: random, QMC, and K-medoids designs on the three-component simplex.
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
n <- 15L # Matches the plotted panels; the paper caption instead says n = 12.
candidate_set <- unique(simplex_trans(spacefillr::generate_sobol_set(10000, p - 1)))
designs <- list(Random = get_rand(candidate_set, n, seed = 3),
                QMC = get_sobol(candidate_set, n, seed = 3),
                `K-medoids` = get_k_medoids(candidate_set, n, seed = 1))
filenames <- c(Random = "rand_X.png", QMC = "sobol_X.png", `K-medoids` = "k_med_X.png")
for (method in names(designs)) {
  save_ternary(filenames[[method]], function() {
    ternary_axes()
    Ternary::TernaryPoints(designs[[method]]$X_sub, pch = 21, bg = "green", lwd = 0.9, cex = 2)
  })
}
saveRDS(list(n = n, candidate_set = candidate_set, designs = designs,
             seeds = c(Random = 3L, QMC = 3L, `K-medoids` = 1L)),
        file.path(output_dir, "Figure_2_results.rds"))
