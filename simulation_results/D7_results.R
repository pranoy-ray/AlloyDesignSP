rm(list = ls())
source("D7_lib.R")

#################################################### D7
D7 = read.csv("RHEA7unique.csv")
X <- as.matrix(D7[,4:10])
X <- scales::rescale(X, to = c(0,1), from = c(min(X),max(X)))
y <- as.matrix(D7[,3], ncol = 1)

n <- c(50,60,70,80,90)
res <- list()
for(i in 1:length(n)) {
  print(i)
  res[[i]] <- simu(X, y, n[i], k = 30)
}
save(res, file = "D7_results.RData")
load("D7_results.RData")

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


plot_test_sets(c(3, 5, 7), "NMAE (%)", "combined_Bulk_Modulus_D7_NMAE.jpeg")
plot_test_sets(c(4, 6, 8), "RMSE", "combined_Bulk_Modulus_D7_RMSE.jpeg")
