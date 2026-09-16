# Figure 6: Runtime Comparison
rm(list = ls())
library(ggplot2)
method_order <- c("Random", "QMC", "K_medoids", "SP", "SSP")
method_colors <- c(
  "Random" = "#D7191C",
  "QMC" = "#1A9641",
  "K_medoids" = "#984EA3",
  "SP" = "#FF7F00",
  "SSP" = "#0072B2"
)
method_linetypes <- c(Random = "solid", QMC = "longdash", K_medoids = "solid",
                      SP = "solid", SSP = "dashed")
method_shapes <- c(Random = 16, QMC = 17, K_medoids = 15, SP = 1, SSP = 4)

#################################################
################################################# D4
n = c(20,30,40,50,60)
load("D4_results.RData")
rand_rt <- c(median(res[[1]]$rt_df$Random), median(res[[2]]$rt_df$Random),
  median(res[[3]]$rt_df$Random), median(res[[4]]$rt_df$Random), median(res[[5]]$rt_df$Random))
sobol_rt <- c(median(res[[1]]$rt_df$Sobol), median(res[[2]]$rt_df$Sobol),
  median(res[[3]]$rt_df$Sobol), median(res[[4]]$rt_df$Sobol), median(res[[5]]$rt_df$Sobol))
k_medoids_rt <- c(median(res[[1]]$rt_df$K_medoids), median(res[[2]]$rt_df$K_medoids),
  median(res[[3]]$rt_df$K_medoids), median(res[[4]]$rt_df$K_medoids), median(res[[5]]$rt_df$K_medoids))
sp_rt <- c(median(res[[1]]$rt_df$SP), median(res[[2]]$rt_df$SP),
  median(res[[3]]$rt_df$SP), median(res[[4]]$rt_df$SP), median(res[[5]]$rt_df$SP))
ssp_rt <- c(median(res[[1]]$rt_df$SSP), median(res[[2]]$rt_df$SSP),
  median(res[[3]]$rt_df$SSP), median(res[[4]]$rt_df$SSP), median(res[[5]]$rt_df$SSP))

rt_df <- data.frame(
  n = rep(n, times = 5),
  runtime = c(
    rand_rt,
    sobol_rt,
    k_medoids_rt,
    sp_rt,
    ssp_rt
  ),
  method = rep(method_order, each = length(n))
)

rt_df$method <- factor(rt_df$method, levels = method_order)

rt_df_D4 <- transform(rt_df, dataset = "D4")


#################################################
################################################# D7
n = c(50,60,70,80,90)
load("D7_results.RData")
rand_rt <- c(median(res[[1]]$rt_df$Random), median(res[[2]]$rt_df$Random),
  median(res[[3]]$rt_df$Random), median(res[[4]]$rt_df$Random), median(res[[5]]$rt_df$Random))
sobol_rt <- c(median(res[[1]]$rt_df$Sobol), median(res[[2]]$rt_df$Sobol),
  median(res[[3]]$rt_df$Sobol), median(res[[4]]$rt_df$Sobol), median(res[[5]]$rt_df$Sobol))
k_medoids_rt <- c(median(res[[1]]$rt_df$K_medoids), median(res[[2]]$rt_df$K_medoids),
  median(res[[3]]$rt_df$K_medoids), median(res[[4]]$rt_df$K_medoids), median(res[[5]]$rt_df$K_medoids))
sp_rt <- c(median(res[[1]]$rt_df$SP), median(res[[2]]$rt_df$SP),
  median(res[[3]]$rt_df$SP), median(res[[4]]$rt_df$SP), median(res[[5]]$rt_df$SP))
ssp_rt <- c(median(res[[1]]$rt_df$SSP), median(res[[2]]$rt_df$SSP),
  median(res[[3]]$rt_df$SSP), median(res[[4]]$rt_df$SSP), median(res[[5]]$rt_df$SSP))

rt_df <- data.frame(
  n = rep(n, times = 5),
  runtime = c(
    rand_rt,
    sobol_rt,
    k_medoids_rt,
    sp_rt,
    ssp_rt
  ),
  method = rep(method_order, each = length(n))
)

rt_df$method <- factor(rt_df$method, levels = method_order)

rt_df <- rbind(rt_df_D4, transform(rt_df, dataset = "D7"))

jpeg("D4_D7_runtime.jpeg", width = 14, height = 6, units = "in", res = 600, quality = 100)
runtime_plot <- ggplot(rt_df, aes(x = n, y = runtime, color = method, group = method)) +
  geom_line(aes(linetype = method), linewidth = 0.85, alpha = 1) +
  geom_point(aes(shape = method), size = 3, alpha = 1) +
  facet_wrap(~dataset, nrow = 1, scales = "free") +
  scale_color_manual(
    values = method_colors,
    breaks = method_order,
    labels = c("Random", "QMC", "K_medoids", "SP", "SSP")
  ) +
  scale_linetype_manual(values = method_linetypes, breaks = method_order) +
  scale_shape_manual(values = method_shapes, breaks = method_order) +
  labs(
    x = "sample size n",
    y = "Runtime (sec)",
    color = "Method", linetype = "Method", shape = "Method"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 22),
    panel.spacing = grid::unit(1, "cm"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(face = "bold", size = 25),
    axis.title.y = element_text(face = "bold", size = 25),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.title = element_text(face = "bold", size = 24),
    legend.text = element_text(size = 22)
  )
print(runtime_plot)
dev.off()
