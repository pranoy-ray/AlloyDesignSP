# Figure 9: D7 coordinate maxima at n = 70.
# Resolve paths for Rscript and source(), including invocation from another directory.
entry_file <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(files)) tail(files, 1)[[1]] else
    sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(entry_file, mustWork = TRUE)), "common.R"), local = TRUE)

require_packages(c(sampling_packages, plot_packages))
library(ggplot2)
source(file.path(script_dir, "sampling_strategy.R"), local = TRUE)

D7 <- read.csv(file.path(script_dir, "RHEA7unique.csv"))
X <- as.matrix(D7[, c("Mo", "Nb", "Ta", "Ti", "V", "W", "Zr")]) / 128
n <- 70
replications <- 30
methods <- list(Random = get_rand, QMC = get_sobol, "K-medoids" = get_k_medoids,
                SP = get_sp, SSP = get_ssp)
colors <- method_colors

# Record the minimum and maximum separately for each coordinate.
records <- list()
design_indices <- list()
for (method in names(methods)) {
  message("Sampling ", method, "...")
  for (i in seq_len(replications)) {
    design <- methods[[method]](X, n, i)
    stopifnot(nrow(design$X_sub) == n, all(is.finite(design$X_sub)),
              length(unique(design$idx)) == n)
    key <- paste(method, i, sep = "_")
    records[[key]] <- data.frame(method = method, replication = i, n = n,
                                coordinate = colnames(design$X_sub),
                                minimum = apply(design$X_sub, 2, min),
                                maximum = apply(design$X_sub, 2, max))
    design_indices[[key]] <- design$idx
    if (i %% 5 == 0) message("  Completed ", i, "/", replications)
  }
}
D7_level_res <- do.call(rbind, records)
rownames(D7_level_res) <- NULL
D7_level_res$method <- factor(D7_level_res$method, levels = names(methods))

save(D7_level_res, design_indices, file = file.path(output_dir, "D7_max_level.Rdata"))
write.csv(D7_level_res, file.path(output_dir, "D7_max_level.csv"), row.names = FALSE)

# Each box summarizes 30 replications; dots show individual designs.
plot_data <- tidyr::pivot_longer(D7_level_res, maximum,
                                 names_to = "statistic", values_to = "value")
plot_data$statistic <- factor(plot_data$statistic, levels = "maximum",
                              labels = "Maximum")
plot_data$coordinate <- factor(plot_data$coordinate, levels = colnames(X))
level_plot <- ggplot(plot_data, aes(method, value, fill = method)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1),
             size = 1.5, alpha = 0.55) +
  facet_grid(statistic ~ coordinate) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(title = "D7: design maxima by coordinate",
       subtitle = "n = 70; 30 replications per method",
       x = NULL, y = "Value") +
  theme_bw(base_size = 22) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        axis.title = element_text(size = 24),
        axis.text.y = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 17),
        strip.text = element_text(face = "bold", size = 22),
        plot.title = element_text(face = "bold", size = 28),
        plot.subtitle = element_text(size = 22))

ggsave(file.path(output_dir, "D7_max_level.jpeg"), plot = level_plot,
       width = 18, height = 5, units = "in", dpi = 600, quality = 100)
