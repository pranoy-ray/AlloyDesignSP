# Regenerate all results and figures using the full paper experiment sizes.
# Resolve paths for Rscript and source(), including invocation from another directory.
entry_file <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(files)) tail(files, 1)[[1]] else
    sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(entry_file, mustWork = TRUE)), "common.R"), local = TRUE)

require_packages(c(sampling_packages, plot_packages, ternary_packages, "rkriging"))
scripts <- c("Figure_2.R", "Figure_3.R", "Figure_4.R", "D4.R", "plot_D4.R",
             "D7.R", "plot_D7.R", "plot_runtime.R", "D7_max_level.R")
for (script in scripts) {
  message("Running ", script)
  source(file.path(script_dir, script), local = new.env(parent = globalenv()))
}
message("All results and figures saved in ", output_dir)
