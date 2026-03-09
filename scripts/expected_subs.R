library(ggplot2)
library(viridis)

sub_rates <- snakemake@params[["sub_rates"]]

max_count <- 20
distribution_df <- sapply(
  X = sub_rates,
  FUN = function(rate) dbinom(x = seq(0, max_count, by = 1), size = 150, prob = rate)
) |>
    as.data.frame()

colnames(distribution_df) <- sub_rates

distribution_df_long <- distribution_df |>
  dplyr::mutate(count = seq(0, max_count, by = 1)) |>
  tidyr::pivot_longer(
    cols = -count,
    names_to = "p_n",
    values_to = "prob_of_count"
  ) |>
  dplyr::mutate(
    count = as.numeric(count),
    p_n = as.factor(p_n)
  )

binom_plot <- ggplot(distribution_df_long, aes(x = count, y = prob_of_count)) +
  geom_line(aes(color = p_n)) +
  scale_color_viridis_d() +
  labs(
    x = "Number of T>C Substitutions",
    y = "Probability of Observing Count",
    color = bquote(p[n])
  ) +
  theme_bw()

ggsave(
  binom_plot,
  filename = snakemake@output[["binom_plot_png"]],
  width = 8,
  height = 6
)
ggsave(
  binom_plot,
  filename = snakemake@output[["binom_plot_pdf"]],
  width = 8,
  height = 6
)