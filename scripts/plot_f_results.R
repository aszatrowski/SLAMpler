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
  aes(x = true_f_new, y = map, fill = read_count)) +
  geom_line(aes(color = read_count)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.25) +
  scale_color_viridis_d(begin = 0.8, end = 0, guide = "none") +
  scale_fill_viridis_d(begin = 0.8, end = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    fill = "Read Count",
    x = bquote("True " ~ p[n]),
    y = bquote("Posterior " ~ p[n]),
  ) +
  theme_bw()
ggsave(
  f_dist_plot,
  filename = snakemake@output[["f_dist_plot_png"]],
  width = snakemake@params[["plot_width"]],
  height = snakemake@params[["plot_height"]]
)
ggsave(
  f_dist_plot,
  filename = snakemake@output[["f_dist_plot_pdf"]],
  width = snakemake@params[["plot_width"]],
  height = snakemake@params[["plot_height"]]
)