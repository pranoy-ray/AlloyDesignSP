# Verify full saved runs without repeating the expensive simulations.
# Resolve paths for Rscript and source(), including invocation from another directory.
entry_file <- local({
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(files)) tail(files, 1)[[1]] else
    sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
})
source(file.path(dirname(normalizePath(entry_file, mustWork = TRUE)), "common.R"), local = TRUE)


source(file.path(script_dir, "MMSE_calculation.R"), local = TRUE)
# Analytical kernel checks: interpolation and the one-training-point solution.
stopifnot(max(MSE(diag(3), diag(3), 0.45)) < 1e-12,
          isTRUE(all.equal(MSE(matrix(c(0, 0), 1), matrix(c(1, 0), 1), 1), 1 - exp(-2))))

expected_sizes <- list(D4 = c(20, 30, 40, 50, 60), D7 = c(50, 60, 70, 80, 90))
all_runtimes <- list()
for (dataset in names(expected_sizes)) {
  saved <- read_results(dataset)
  stopifnot(identical(saved$metadata$sizes, expected_sizes[[dataset]]),
            saved$metadata$replications == 30L, saved$metadata$kappa == 1,
            identical(saved$metadata$seeds, 1:30))
  data <- read.csv(file.path(script_dir, if (dataset == "D4") "pca_pspall4.csv" else "RHEA7unique.csv"))
  X <- as.matrix(data[, if (dataset == "D4") c("Al", "Nb", "Ti", "Zr") else
                       c("Mo", "Nb", "Ta", "Ti", "V", "W", "Zr")]) / 128
  stopifnot(all(abs(rowSums(X) - 1) < 1e-12), all(X >= 0 & X <= 1))
  rows <- list()
  suffixes <- if (dataset == "D4") c("1", "2") else ""
  for (method in names(method_objects)) {
    results <- saved[[paste(dataset, method_objects[[method]], sep = "_")]]
    indices <- saved$design_indices[[method]]
    stopifnot(length(results) == 5L, all(lengths(results) == 30L),
              length(indices) == 5L, all(lengths(indices) == 30L))
    for (i in seq_along(results)) for (j in seq_len(30L)) {
      metrics <- results[[i]][[j]]
      train <- indices[[i]][[j]]
      stopifnot(length(train) == expected_sizes[[dataset]][i], !anyDuplicated(train),
                all(train %in% seq_len(nrow(X))), all(is.finite(unlist(metrics))),
                all(unlist(metrics) >= 0))
      test <- setdiff(seq_len(nrow(X)), train)
      boundary <- rowSums(X[test, , drop = FALSE] == 0 | X[test, , drop = FALSE] == 1) > 0
      weight <- mean(boundary)
      stopifnot(weight > 0, weight < 1)
      for (suffix in suffixes) {
        # Full-set errors must be the size-weighted combination of its disjoint subsets.
        get_metric <- function(subset, metric) metrics[[paste0(subset, "_", metric, suffix)]]
        stopifnot(isTRUE(all.equal(get_metric("full", "mae"),
          weight * get_metric("bound", "mae") + (1 - weight) * get_metric("inter", "mae"), tolerance = 1e-10)),
          isTRUE(all.equal(get_metric("full", "rmse")^2,
          weight * get_metric("bound", "rmse")^2 + (1 - weight) * get_metric("inter", "rmse")^2, tolerance = 1e-10)))
      }
      rows[[length(rows) + 1L]] <- data.frame(method = method, n = expected_sizes[[dataset]][i],
                                             replication = j, as.data.frame(metrics))
    }
  }
  expected <- do.call(rbind, rows)
  csv <- read.csv(file.path(output_dir, paste0(dataset, "_metrics.csv")))
  stopifnot(nrow(csv) == 750L, isTRUE(all.equal(csv, expected, tolerance = 1e-12)))
  all_runtimes[[dataset]] <- data.frame(dataset = dataset, expected[, c("method", "n", "runtime")])
  message(dataset, ": 750 designs and their full/boundary/interior metrics verified")
}

actual_runtime <- read.csv(file.path(output_dir, "D4_D7_runtime.csv"))
expected_runtime <- aggregate(runtime ~ dataset + method + n, do.call(rbind, all_runtimes), median)
sort_runtime <- function(x) x[order(x$dataset, x$method, x$n), ]
a <- sort_runtime(actual_runtime)
b <- sort_runtime(expected_runtime)
rownames(a) <- rownames(b) <- NULL
stopifnot(isTRUE(all.equal(a, b, tolerance = 1e-12)))

f2 <- readRDS(file.path(output_dir, "Figure_2_results.rds"))
f3 <- readRDS(file.path(output_dir, "Figure_3_results.rds"))
f4 <- readRDS(file.path(output_dir, "Figure_4_results.rds"))
stopifnot(f2$n == 15L, f3$n == 15L, f4$n == 15L, f3$theta == 0.45, f3$kappa == 1)
# The MSE test set contains 40,000 simplex points followed by all three vertices.
stopifnot(identical(dim(f3$test_set), c(40003L, 3L)),
          all(is.finite(f3$test_set)), all(f3$test_set >= 0 & f3$test_set <= 1),
          max(abs(rowSums(f3$test_set) - 1)) < 1e-12,
          identical(unname(tail(f3$test_set, 3L)), diag(3L)))
for (design in c(f2$designs, f3$designs, list(f4$SP, f4$SSP))) {
  stopifnot(nrow(design$X_sub) == 15L, !anyDuplicated(design$idx),
            max(abs(rowSums(design$X_sub) - 1)) < 1e-12)
}
for (method in names(f3$designs)) {
  values <- MSE(f3$designs[[method]]$X_sub, f3$test_set, f3$theta)
  stopifnot(isTRUE(all.equal(values, f3$mse[[method]], tolerance = 1e-12)),
            all(values >= 0 & values <= 1), MMSE(f3$designs[[method]]$X_sub, f3$test_set, f3$theta) == max(values))
}
stopifnot(identical(f3$designs$SP$idx, f4$SP$idx), identical(f3$designs$SSP$idx, f4$SSP$idx))

levels <- new.env()
load(file.path(output_dir, "D7_max_level.Rdata"), envir = levels)
level_csv <- read.csv(file.path(output_dir, "D7_max_level.csv"))
stopifnot(nrow(level_csv) == 1050L, length(levels$design_indices) == 150L,
          all(level_csv$minimum == 0), all(level_csv$maximum <= 1))
# Figure 9 must describe the same seeded n=70 designs used for Figure 8.
saved_d7 <- read_results("D7")
for (method in names(method_objects)) for (j in seq_len(30L)) {
  idx <- levels$design_indices[[paste(method, j, sep = "_")]]
  stopifnot(identical(idx, saved_d7$design_indices[[method]][[3]][[j]]))
  rows <- level_csv[level_csv$method == method & level_csv$replication == j, ]
  design <- X[idx, , drop = FALSE] # X from the D7 loop above.
  stopifnot(identical(rows$coordinate, colnames(design)),
            isTRUE(all.equal(rows$maximum, unname(apply(design, 2, max)))),
            isTRUE(all.equal(rows$minimum, unname(apply(design, 2, min)))))
}

images <- c("rand_X.png", "sobol_X.png", "k_med_X.png", "MSE_ternary_SP.png",
            "MSE_ternary_SSP_kappa_1.png", "all_X.png", "D4_D7_runtime.jpeg", "D7_max_level.jpeg",
            outer(c("combined_Bulk_Modulus_D4", "combined_Formation_Energy_D4", "combined_Bulk_Modulus_D7"),
                  c("_NMAE.jpeg", "_RMSE.jpeg"), paste0))
require_packages("magick")
for (filename in images) {
  info <- magick::image_info(magick::image_read(file.path(output_dir, filename)))
  stopifnot(nrow(info) == 1, info$width > 1000, info$height > 1000)
}
message("Verified 1,500 prediction designs, 150 coordinate-maxima designs, and all 14 images.")
