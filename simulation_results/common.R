# Shared paths and dependency checks. Entry scripts resolve this file relative
# to their own location, so neither setwd() nor a particular shell directory is needed.
common_file <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(files)) tail(files, 1)[[1]] else
    sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
})
script_dir <- dirname(normalizePath(common_file, mustWork = TRUE))
output_dir <- file.path(script_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         ". Install them with install.packages(c(",
         paste(sprintf('"%s"', missing), collapse = ", "), ")).", call. = FALSE)
  }
}

sampling_packages <- c("spacefillr", "SPlit", "cluster", "support")
plot_packages <- c("ggplot2", "tidyr", "tidyselect")
ternary_packages <- c("Ternary", "magick", "PlotTools")

method_colors <- c(Random = "#D7191C", QMC = "#1A9641", "K-medoids" = "#984EA3",
                   SP = "#FF7F00", SSP = "#0072B2")
method_objects <- c(Random = "rand_res", QMC = "qmc_res", "K-medoids" = "k_res",
                    SP = "sp_res", SSP = "ssp_res")

read_results <- function(dataset) {
  path <- file.path(output_dir, paste0(dataset, ".Rdata"))
  if (!file.exists(path)) stop("Missing ", path, "; run ", dataset, ".R first.")
  saved <- new.env(parent = emptyenv())
  load(path, envir = saved)
  if (!isTRUE(saved$metadata$complete)) {
    stop(path, " is incomplete; finish running ", dataset, ".R before plotting.")
  }
  saved
}

ternary_axes <- function() {
  par(mar = c(0.6, 0.6, 0.6, 4.2), xpd = NA)
  Ternary::TernaryPlot(
    atip = "x1", btip = "x2", ctip = "x3", lab.cex = 1.7, axis.cex = 1.5,
    grid.lines = 5, grid.minor.lines = 1, grid.lty = "solid",
    grid.minor.lty = "dotted", col = NA, grid.col = "grey85",
    axis.col = "grey40", ticks.col = "grey40", axis.labels = seq(0, 1, by = 0.2),
    axis.rotate = FALSE, padding = 0.08
  )
}

save_ternary <- function(filename, draw) {
  path <- file.path(output_dir, filename)
  png(path, width = 7, height = 6, units = "in", res = 600)
  tryCatch(draw(), finally = dev.off())
  magick::image_write(magick::image_trim(magick::image_read(path)), path = path)
  invisible(path)
}
