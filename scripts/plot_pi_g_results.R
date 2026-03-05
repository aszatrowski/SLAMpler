library(dplyr)
library(ggplot2)

pi_g_samples <- readr::read_csv(unlist(snakemake@input[["pi_samples"]]))

calibration_df <- pi_g_samples |>
  group_by(true_pi) |>
  summarise(
    mean_estimated_pi = mean(pi_g),
    lower_ci = quantile(pi_g, 0.025),
    upper_ci = quantile(pi_g, 0.975)
  )

pi_g_dist_plot <- ggplot(calibration_df, aes(x = true_pi, y = mean_estimated_pi, color = true_pi)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.005) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#434343") +
  scale_color_viridis_c(begin = 0, end = 0.6) +
  labs(
    color = bquote(pi[g]),
    x = bquote("True " ~ pi[g]),
    y = bquote("Posterior " ~ pi[g])
  ) +
  theme_bw()
ggsave(
  pi_g_dist_plot,
  filename = snakemake@output[["pi_dist_plot"]],
  width = 8,
  height = 6
)