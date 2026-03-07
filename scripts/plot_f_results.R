library(dplyr)
library(ggplot2)

f_samples <- readr::read_csv(unlist(snakemake@input[["pi_samples"]]))

calibration_df <- f_samples |>
  mutate(
    read_count = as.factor(read_count)
  ) |>
  group_by(true_f_new, read_count) |>
  summarise(
    # mean a posteriori estimate
    map = mean(f_new),
    # credible interval
    lower_ci = quantile(f_new, 0.025),
    upper_ci = quantile(f_new, 0.975)
  )

f_dist_plot <- ggplot(
  calibration_df,
  aes(x = true_f_new, y = map, color = read_count, group = read_count)
) +
  geom_point(
    position = position_dodge(width = 0.005)
  ) +
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci),
    width = 0.005,
    position = position_dodge(width = 0.005)
  ) +
  scale_color_viridis_d(begin = 0, end = 0.8) +
  labs(
    color = "Read Coverage",
    x = bquote("True " ~ p[n]),
    y = bquote("Posterior " ~ p[n]),
  ) +
  theme_bw()
ggsave(
  f_dist_plot,
  filename = snakemake@output[["f_dist_plot"]],
  width = 8,
  height = 6
)