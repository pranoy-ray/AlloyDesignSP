# Shared D4/D7 evaluation. Sampling functions are loaded by the entry script.
predictive_metrics <- function(X, responses, train_idx) {
  test_idx <- setdiff(seq_len(nrow(X)), train_idx)
  boundary <- rowSums(X[test_idx, , drop = FALSE] == 0 |
                      X[test_idx, , drop = FALSE] == 1) > 0
  subsets <- list(full = seq_along(test_idx), bound = which(boundary),
                  inter = which(!boundary))
  stopifnot(all(lengths(subsets) > 0L))
  metrics <- list()
  for (response_id in seq_along(responses)) {
    response <- if (length(responses) == 1L) "" else names(responses)[response_id]
    y <- responses[[response_id]]
    fit <- rkriging::Fit.Kriging(X[train_idx, , drop = FALSE], y[train_idx],
                                kernel.parameters = list(type = "Gaussian"))
    prediction <- rkriging::Predict.Kriging(fit, X[test_idx, , drop = FALSE])$mean
    errors <- y[test_idx] - as.vector(prediction)
    stopifnot(all(is.finite(errors)), mean(abs(y)) > 0)
    for (test_set in names(subsets)) {
      e <- errors[subsets[[test_set]]]
      # Paper Eq. (18): normalize by the entire candidate set; report percent.
      metrics[[paste0(test_set, "_mae", response)]] <- 100 * mean(abs(e)) / mean(abs(y))
      metrics[[paste0(test_set, "_rmse", response)]] <- sqrt(mean(e^2))
    }
  }
  metrics
}

run_experiment <- function(dataset, X, responses, sizes, replications = 30L) {
  stopifnot(all(is.finite(X)), all(abs(rowSums(X) - 1) < 1e-10),
            all(lengths(responses) == nrow(X)), all(sizes < nrow(X)))
  methods <- list(Random = get_rand, QMC = get_sobol, "K-medoids" = get_k_medoids,
                  SP = get_sp, SSP = get_ssp)
  saved <- new.env(parent = emptyenv())
  saved$metadata <- list(dataset = dataset, sizes = sizes, replications = replications,
                         seeds = seq_len(replications), kappa = 1, complete = FALSE,
                         started = Sys.time(), session = sessionInfo())
  saved$design_indices <- setNames(vector("list", length(methods)), names(methods))
  rows <- list()
  data_file <- file.path(output_dir, paste0(dataset, ".Rdata"))
  for (method in names(methods)) {
    object <- paste(dataset, method_objects[[method]], sep = "_")
    results <- vector("list", length(sizes))
    indices <- vector("list", length(sizes))
    for (i in seq_along(sizes)) {
      results[[i]] <- vector("list", replications)
      indices[[i]] <- vector("list", replications)
      for (j in seq_len(replications)) {
        runtime <- system.time(design <- methods[[method]](X, sizes[i], j))[["elapsed"]]
        stopifnot(nrow(design$X_sub) == sizes[i],
                  length(unique(design$idx)) == sizes[i],
                  all(is.finite(design$X_sub)))
        ans <- predictive_metrics(X, responses, design$idx)
        ans$runtime <- runtime
        results[[i]][[j]] <- ans
        indices[[i]][[j]] <- design$idx
        rows[[length(rows) + 1L]] <- data.frame(method = method, n = sizes[i],
                                               replication = j, as.data.frame(ans))
        if (j %% 5L == 0L) message(dataset, " / ", method, " / n = ", sizes[i],
                                    ": ", j, "/", replications)
      }
      assign(object, results, envir = saved)
      saved$design_indices[[method]] <- indices
      # Keep a checkpoint after each design size; plotters reject incomplete runs.
      save(list = ls(saved), envir = saved, file = data_file)
    }
  }
  saved$metadata$complete <- TRUE
  saved$metadata$finished <- Sys.time()
  save(list = ls(saved), envir = saved, file = data_file)
  write.csv(do.call(rbind, rows), file.path(output_dir, paste0(dataset, "_metrics.csv")),
            row.names = FALSE)
  writeLines(trimws(capture.output(sessionInfo()), which = "right"),
             file.path(output_dir, paste0(dataset, "_session.txt")))
  message("Saved ", data_file)
  invisible(saved)
}
