# Prediction error plots for D7.
# Resolve paths for Rscript and source(), including invocation from another directory.
entry_file <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(files)) tail(files, 1)[[1]] else
    sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(entry_file, mustWork = TRUE)), "common.R"), local = TRUE)

require_packages(plot_packages)
library(ggplot2)
saved <- read_results("D7")
sizes <- saved$metadata$sizes
replications <- saved$metadata$replications
objects <- c(Random = "D7_rand_res", QMC = "D7_qmc_res",
             "K-medoids" = "D7_k_res", SP = "D7_sp_res", SSP = "D7_ssp_res")
colors <- method_colors
metric_names <- as.vector(outer(c("full_", "bound_", "inter_"),
                               c("mae", "rmse"), paste0))

# One row per method, design size, and replication; values are used unchanged.
rows <- lapply(names(objects), function(method) {
  results <- get(objects[[method]], envir = saved, inherits = FALSE)
  stopifnot(length(results) == length(sizes), all(lengths(results) == replications))
  do.call(rbind, lapply(seq_along(sizes), function(i) {
    values <- do.call(rbind, lapply(results[[i]], as.data.frame))
    stopifnot(setequal(names(values), c(metric_names, "runtime")),
              all(is.finite(as.matrix(values))))
    values <- values[, metric_names, drop = FALSE]
    data.frame(method = method, n = sizes[i], replication = seq_len(replications), values)
  }))
})
plot_data <- tidyr::pivot_longer(do.call(rbind, rows), cols = tidyselect::all_of(metric_names),
  names_to = c("test_set", "metric"),
  names_pattern = "(full|bound|inter)_(mae|rmse)", values_to = "value")
plot_data$method <- factor(plot_data$method, levels = names(objects))
plot_data$n <- factor(plot_data$n, levels = sizes)
plot_data$test_set <- factor(plot_data$test_set, levels = c("full", "bound", "inter"),
                             labels = c("Full testing set", "Boundary only", "Interior only"))

save_boxplot <- function(metric_id) {
  data <- plot_data[plot_data$metric == metric_id, ]
  metric_label <- if (metric_id == "mae") "NMAE" else "RMSE"
  p <- ggplot(data, aes(n, value, fill = method)) +
    geom_boxplot(width = 0.75, position = position_dodge(width = 0.85),
                 outlier.shape = 1, outlier.size = 1.2) +
    facet_wrap(~test_set, nrow = 1) +
    scale_fill_manual(values = colors) +
    labs(title = "D7: Bulk Modulus",
         x = "Sample Size n", y = if (metric_id == "mae") "NMAE (%)" else "RMSE (GPa)", fill = "Method") +
    guides(fill = guide_legend(nrow = 1)) +
    theme_bw(base_size = 16) +
    theme(axis.title = element_text(face = "bold", size = 21),
          axis.text = element_text(size = 15), strip.text = element_text(face = "bold", size = 19),
          legend.position = "bottom", legend.title = element_text(face = "bold"),
          legend.text = element_text(size = 14), plot.title = element_text(face = "bold"))
  filename <- paste0("combined_Bulk_Modulus_D7_", metric_label, ".jpeg")
  ggsave(file.path(output_dir, filename), plot = p, width = 14, height = 7,
         units = "in", dpi = 600, quality = 100)
  invisible(p)
}

save_boxplot("mae")
save_boxplot("rmse")
