library(dplyr)
library(ggplot2)

sub_rate_new <- as.numeric(snakemake@wildcards[["sub_rate_new"]])
pi_g_samples <- readr::read_csv(unlist(snakemake@input[["pi_samples"]]), show_col_types = FALSE)

calibration_df <- pi_g_samples |>
  mutate(
    read_count = as.factor(read_count)
  ) |>
  group_by(true_pi, read_count) |>
  summarise(
    map = mean(pi_g),
    lower_ci = quantile(pi_g, 0.025),
    upper_ci = quantile(pi_g, 0.975),
    .groups = "drop_last"
  )

pi_g_dist_plot <- ggplot(
  calibration_df,
  aes(x = true_pi, y = map, fill = read_count)) +
  geom_line(aes(color = read_count)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.15) +
  scale_color_viridis_d(begin = 0.8, end = 0, guide = "none") +
  scale_fill_viridis_d(begin = 0.8, end = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    fill = "Read Count",
    x = bquote("True " ~ pi[g]),
    y = bquote("Posterior " ~ pi[g]),
    caption = bquote("Substitution Rate:" ~ .(sub_rate_new))
  ) +
  theme_bw()

ggsave(
  pi_g_dist_plot,
  filename = snakemake@output[["pi_dist_plot_png"]],
  width = 10,
  height = 6
)
ggsave(
  pi_g_dist_plot,
  filename = snakemake@output[["pi_dist_plot_pdf"]],
  width = 10,
  height = 6
)