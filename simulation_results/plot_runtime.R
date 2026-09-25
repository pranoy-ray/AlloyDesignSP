# Figure 5: median sampling-only runtimes.
# Resolve paths for Rscript and source(), including invocation from another directory.
entry_file <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(files)) tail(files, 1)[[1]] else
    sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(entry_file, mustWork = TRUE)), "common.R"), local = TRUE)

require_packages("ggplot2")
library(ggplot2)
objects <- method_objects
colors <- method_colors
saved_results <- list(D4 = read_results("D4"), D7 = read_results("D7"))
sizes <- lapply(saved_results, function(x) x$metadata$sizes)
linetypes <- c(Random = "solid", QMC = "longdash", "K-medoids" = "solid",
               SP = "solid", SSP = "dashed")
shapes <- c(Random = 16, QMC = 17, "K-medoids" = 15, SP = 1, SSP = 4)

# Read the sampling-only elapsed seconds recorded in each replication.
rows <- lapply(names(sizes), function(dataset) {
  saved <- saved_results[[dataset]]
  replications <- saved$metadata$replications
  do.call(rbind, lapply(names(objects), function(method) {
    results <- get(paste(dataset, objects[[method]], sep = "_"),
                   envir = saved, inherits = FALSE)
    stopifnot(length(results) == length(sizes[[dataset]]), all(lengths(results) == replications))
    do.call(rbind, lapply(seq_along(sizes[[dataset]]), function(i) {
      runtime <- unlist(lapply(results[[i]], `[[`, "runtime"), use.names = FALSE)
      stopifnot(length(runtime) == replications, all(is.finite(runtime)), all(runtime >= 0))
      data.frame(dataset = dataset, method = method, n = sizes[[dataset]][i],
                 replication = seq_len(replications), runtime = runtime)
    }))
  }))
})
runtime_data <- do.call(rbind, rows)
runtime_data$method <- factor(runtime_data$method, levels = names(objects))

# Paper Figure 5: each point is the median of the 30 replications.
runtime_summary <- aggregate(runtime ~ dataset + method + n, runtime_data, median)
write.csv(runtime_summary, file.path(output_dir, "D4_D7_runtime.csv"), row.names = FALSE)
runtime_plot <- ggplot(runtime_summary, aes(n, runtime, color = method, group = method)) +
  geom_line(aes(linetype = method), linewidth = 0.85) +
  geom_point(aes(shape = method, size = method)) +
  facet_wrap(~dataset, nrow = 1, scales = "free") +
  scale_x_continuous(breaks = sort(unique(unlist(sizes)))) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  scale_shape_manual(values = shapes) +
  scale_size_manual(name = "Method", values = c(Random = 8, QMC = 7,
                    "K-medoids" = 7, SP = 7, SSP = 7)) +
  labs(x = "Sample Size n", y = "Median runtime (sec)",
       color = "Method", linetype = "Method", shape = "Method") +
  theme_bw() +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 22),
        panel.spacing = grid::unit(1, "cm"), panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold", size = 25), axis.text = element_text(size = 20),
        legend.title = element_text(face = "bold", size = 24), legend.text = element_text(size = 22))

ggsave(file.path(output_dir, "D4_D7_runtime.jpeg"), plot = runtime_plot,
       width = 14, height = 6, units = "in", dpi = 600, quality = 100)
