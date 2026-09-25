# D7: 30 replications per method and design size.
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

data <- read.csv(file.path(script_dir, "RHEA7unique.csv"))
# Both supplied datasets contain atom counts summing to 128 per composition.
X <- as.matrix(data[, c("Mo", "Nb", "Ta", "Ti", "V", "W", "Zr")]) / 128
run_experiment("D7", X, list(Bulk = data$Bulk), sizes = c(50, 60, 70, 80, 90))
