library(dplyr)
library(ggplot2)

pi_g_samples <- readr::read_csv(unlist(snakemake@input[["pi_samples"]]))

calibration_df <- pi_g_samples |>
  mutate(
    read_count = as.factor(read_count)
  ) |>
  group_by(true_pi, read_count) |>
  summarise(
    # mean a posteriori estimate
    map = mean(pi_g),
    # credible interval
    lower_ci = quantile(pi_g, 0.025),
    upper_ci = quantile(pi_g, 0.975)
  )

pi_g_dist_plot <- ggplot(
  calibration_df,
  aes(x = true_pi, y = map, color = read_count, group = read_count)
) +
  geom_point(
    position = position_dodge(width = 0.0125)
  ) +
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci),
    width = 0.005,
    position = position_dodge(width = 0.0125)
  ) +
  scale_color_viridis_d(begin = 0, end = 0.8) +
  labs(
    color = "reads",
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