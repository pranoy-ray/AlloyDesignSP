# D4: 30 replications per method and design size.
# Resolve paths for Rscript and source(), including invocation from another directory.
entry_file <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(files)) tail(files, 1)[[1]] else
    sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(entry_file, mustWork = TRUE)), "common.R"), local = TRUE)

require_packages(c(sampling_packages, "rkriging"))
source(file.path(script_dir, "sampling_strategy.R"), local = TRUE)
source(file.path(script_dir, "experiment_helpers.R"), local = TRUE)

data <- read.csv(file.path(script_dir, "pca_pspall4.csv"))
# Both supplied datasets contain atom counts summing to 128 per composition.
X <- as.matrix(data[, c("Al", "Nb", "Ti", "Zr")]) / 128
run_experiment("D4", X, list(`1` = data$formation_energy, `2` = data$bulk_modulus), sizes = c(20, 30, 40, 50, 60))
