rm(list = ls())
source("D4_lib.R")

#################################################### D4
D4 = read.csv("pca_pspall4.csv")
X <- as.matrix(D4[,c(54,55,56,57)])
X <- scales::rescale(X, to = c(0,1), from = c(min(X),max(X)))
y <- as.matrix(D4[,c(52,53)])

n <- c(20,30,40,50,60)
res <- list()
for(i in 1:length(n)) {
  print(i)
  res[[i]] <- simu(X, y, n[i], k = 30)
}
save(res, file = "D4_results.RData")
load("D4_results.RData")

# duplication check
purrr::map(res, "duplicate_df") # no sampling method produces duplicated points

# Compare full, boundary, and interior testing sets with one shared legend.
plot_test_sets <- function(property_ids, ylab, file_name) {
  test_sets <- c("Full testing set", "Boundary only", "Interior only")
  plots <- lapply(property_ids, function(id) {
    plot_property_boxplot(res, n, property_id = id, ylab = ylab,
                          method_order = c("Random", "QMC", "K_medoids", "SP", "SSP"))
  })
  plot_df <- dplyr::bind_rows(lapply(seq_along(plots), function(i) {
    transform(plots[[i]]$data, test_set = factor(test_sets[i], levels = test_sets))
  }))
  combined_plot <- plots[[1]] + plot_df +
    facet_wrap(~test_set, nrow = 1) +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold", size = 24),
          legend.text = element_text(size = 22),
          strip.text = element_text(face = "bold", size = 22),
          panel.spacing = grid::unit(0.25, "cm"))
  jpeg(file_name, width = 10.5, height = 6, units = "in", res = 600, quality = 100)
  on.exit(dev.off())
  print(combined_plot)
}

plot_test_sets(c(3, 7, 11), "NMAE (%)", "combined_Formation_Energy_D4_NMAE.jpeg")
plot_test_sets(c(4, 8, 12), "RMSE", "combined_Formation_Energy_D4_RMSE.jpeg")
plot_test_sets(c(5, 9, 13), "NMAE (%)", "combined_Bulk_Modulus_D4_NMAE.jpeg")
plot_test_sets(c(6, 10, 14), "RMSE", "combined_Bulk_Modulus_D4_RMSE.jpeg")
