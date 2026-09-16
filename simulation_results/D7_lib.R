# Library for Section 4: Results and Discussion. D7
library(ggplot2)
library(dplyr)
library(tidyr)
source("sampling_strategy.R")

# Generate training and testing set using different sampling strategies.
train_test <- function(X, method, n, seed) {
  duplicate = FALSE
  if (method == "sp") {
    rt <- system.time({idx <- get_sp_idx(X, n, seed)})[3]
  } else if (method == "sobol") {
    rt <- system.time({idx <- get_sobol_idx(X, n, seed)})[3]
  } else if (method == "rand") {
    rt <- system.time({idx <- get_rand_idx(X, n, seed)})[3]
  } else if (method == "ssp") {
    rt <- system.time({idx <- get_ssp_idx(X, n, seed)})[3]
  } else if (method == "K_medoids") {
    rt <- system.time({idx <- get_k_medoids_idx(X, n, seed)})[3]
  }
  if(length(unique(idx)) != n) {duplicate = TRUE}

  idx_2 <- setdiff(1:nrow(X), idx)
  boundary_idx <- which(rowSums(X == 0 | X == 1) > 0)
  test_boundary_idx = intersect(idx_2, boundary_idx)
  test_interior_idx = setdiff(idx_2, test_boundary_idx)

  ans <- list(train_idx = idx, test_idx = idx_2, 
    test_boundary_idx = test_boundary_idx, test_interior_idx = test_interior_idx, runtime = rt, duplicate = duplicate)
  return(ans)
}

# predictive performance
prediction_res <- function(train_X, train_y, test_X, test_y, y) {
  NMAE <- c()
  RMSE <- c()
  fit <- rkriging::Fit.Kriging(train_X, train_y, kernel.parameters = list(type = "Gaussian"))
  pred <- rkriging::Predict.Kriging(fit, test_X)$mean
  mae <- mean(abs(test_y - pred))
  NMAE <- c(NMAE, (100*mae) / mean(abs(y)))
  RMSE <- c(RMSE, sqrt(mean((test_y - pred)^2)))
  ans <- list(NMAE = NMAE, RMSE = RMSE)
  return(ans)
}


# Simulation results across k replications for each sampling method and two response variables: Formation Energy (y1) and Bulk Modulus (y2).
simu_each <- function(X, y, n, k = 30, method) {
  all_y1_nmae <- all_y1_rmse <- c()
  boundary_y1_nmae <- boundary_y1_rmse <- c()
  interior_y1_nmae <- interior_y1_rmse <- c()
  rt <- c()
  duplicate <- c()
  fill_dist <- c()
  
  for(i in 1:k) {
    fill_dist <- c(fill_dist, OSFD::mMdist(X[a$train_idx,], X))
    a <- train_test(X, method, n, i)
    rt <- c(rt, a$runtime)
    duplicate <- c(duplicate, a$duplicate)
    all_res <- prediction_res(X[a$train_idx,], y[a$train_idx,], X[a$test_idx,], y[a$test_idx,], y)
    boundary_res <- prediction_res(X[a$train_idx,], y[a$train_idx,], X[a$test_boundary_idx,], y[a$test_boundary_idx,], y)
    interior_res <- prediction_res(X[a$train_idx,], y[a$train_idx,], X[a$test_interior_idx,], y[a$test_interior_idx,], y)
    all_y1_nmae <- c(all_y1_nmae, all_res$NMAE[1])
    all_y1_rmse <- c(all_y1_rmse, all_res$RMSE[1])
    boundary_y1_nmae <- c(boundary_y1_nmae, boundary_res$NMAE[1])
    boundary_y1_rmse <- c(boundary_y1_rmse, boundary_res$RMSE[1])
    interior_y1_nmae <- c(interior_y1_nmae, interior_res$NMAE[1])
    interior_y1_rmse <- c(interior_y1_rmse, interior_res$RMSE[1])
  }
  ans <- list(
    all_y1_nmae=all_y1_nmae, all_y1_rmse=all_y1_rmse,
    boundary_y1_nmae=boundary_y1_nmae, boundary_y1_rmse=boundary_y1_rmse,
    interior_y1_nmae=interior_y1_nmae, interior_y1_rmse=interior_y1_rmse,
    runtime = rt, duplicate = duplicate)
  return(ans)
}


# simulation: two responses (y1: Formation Energy, y2: Bulk Modulus)
simu <- function(X, y, n, k = 30) {
  sp_res <- simu_each(X, y, n, k, "sp")
  sobol_res <- simu_each(X, y, n, k, "sobol")
  rand_res <- simu_each(X, y, n, k, "rand")
  ssp_res <- simu_each(X, y, n, k, "ssp")
  K_medoids_res <- simu_each(X, y, n, k, "K_medoids")

  rt_df <- data.frame(Random = rand_res$runtime, Sobol = sobol_res$runtime, K_medoids = K_medoids_res$runtime, SP = sp_res$runtime, SSP = ssp_res$runtime)
  duplicate_df <- data.frame(Random = rand_res$duplicate, Sobol = sobol_res$duplicate, K_medoids = K_medoids_res$duplicate, SP = sp_res$duplicate, SSP = ssp_res$duplicate)
  all_y1_nmae_df <- data.frame(Random = rand_res$all_y1_nmae, Sobol = sobol_res$all_y1_nmae, K_medoids = K_medoids_res$all_y1_nmae, SP = sp_res$all_y1_nmae, SSP = ssp_res$all_y1_nmae)
  all_y1_rmse_df <- data.frame(Random = rand_res$all_y1_rmse, Sobol = sobol_res$all_y1_rmse, K_medoids = K_medoids_res$all_y1_rmse, SP = sp_res$all_y1_rmse, SSP = ssp_res$all_y1_rmse)
  boundary_y1_nmae_df <- data.frame(Random = rand_res$boundary_y1_nmae, Sobol = sobol_res$boundary_y1_nmae, K_medoids = K_medoids_res$boundary_y1_nmae, SP = sp_res$boundary_y1_nmae, SSP = ssp_res$boundary_y1_nmae)
  boundary_y1_rmse_df <- data.frame(Random = rand_res$boundary_y1_rmse, Sobol = sobol_res$boundary_y1_rmse, K_medoids = K_medoids_res$boundary_y1_rmse, SP = sp_res$boundary_y1_rmse, SSP = ssp_res$boundary_y1_rmse)
  interior_y1_nmae_df <- data.frame(Random = rand_res$interior_y1_nmae, Sobol = sobol_res$interior_y1_nmae, K_medoids = K_medoids_res$interior_y1_nmae, SP = sp_res$interior_y1_nmae, SSP = ssp_res$interior_y1_nmae)
  interior_y1_rmse_df <- data.frame(Random = rand_res$interior_y1_rmse, Sobol = sobol_res$interior_y1_rmse, K_medoids = K_medoids_res$interior_y1_rmse, SP = sp_res$interior_y1_rmse, SSP = ssp_res$interior_y1_rmse)
  
  ans <- list(
    rt_df=rt_df, duplicate_df=duplicate_df, 
    all_y1_nmae_df=all_y1_nmae_df, all_y1_rmse_df=all_y1_rmse_df,
    boundary_y1_nmae_df=boundary_y1_nmae_df, boundary_y1_rmse_df=boundary_y1_rmse_df,
    interior_y1_nmae_df=interior_y1_nmae_df, interior_y1_rmse_df=interior_y1_rmse_df)
  return(ans)
}


plot_property_boxplot <- function(
  res,
  n,
  property_id,
  ylab = "",
  method_order = NULL,
  methods_keep = NULL,
  ylim = NULL
) {
  plot_df <- lapply(seq_along(n), function(i) {
    df <- res[[i]][[property_id]]
    df$n <- n[i]
    df
  }) |>
    dplyr::bind_rows() |>
    tidyr::pivot_longer(
      cols = -n,
      names_to = "method",
      values_to = "value"
    ) |>
    dplyr::mutate(
      method = dplyr::recode(method, Sobol = "QMC")
    )

  if (!is.null(methods_keep)) {
    methods_keep <- dplyr::recode(methods_keep, Sobol = "QMC")

    plot_df <- plot_df |>
      dplyr::filter(method %in% methods_keep)
  }

  if (!is.null(method_order)) {
    method_order <- dplyr::recode(method_order, Sobol = "QMC")
    plot_df$method <- factor(plot_df$method, levels = method_order)
  }

  method_colors <- c(
    "Random" = "#D7191C",
    "QMC" = "#1A9641",
    "K_medoids" = "#984EA3",
    "SP" = "#FF7F00",
    "SSP" = "#A6CEE3"
  )

  ggplot2::ggplot(plot_df, ggplot2::aes(x = factor(n), y = value, fill = method)) +
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.7,
      outlier.shape = 1
    ) +
    ggplot2::scale_fill_manual(values = method_colors) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(
      x = "Sample Size n",
      y = ylab,
      fill = "Method"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(face = "bold", size = 25),
      axis.title.y = ggplot2::element_text(face = "bold", size = 25),
      axis.text.x = ggplot2::element_text(size = 20),
      axis.text.y = ggplot2::element_text(size = 20),
      legend.title = ggplot2::element_text(face = "bold", size = 20),
      legend.text = ggplot2::element_text(size = 18)
    )
}
