library(dplyr)
library(ggplot2)
library(viridis)

z_hist_tbl <- readr::read_csv(snakemake@input[[1]], show_col_types = FALSE)

z_hist_tbl <- z_hist_tbl |>
  mutate(
    read_id_num = as.numeric(gsub("z\\.", "", read_id)),
    class = case_when(
        class == "old" ~ "Old RNA",
        class == "new" ~ "New RNA",
        TRUE ~ as.character(class)
    )
  )

z_plot <- ggplot(z_hist_tbl, aes(x = iter, y = read_id_num, fill = class)) +
  geom_tile(alpha = 0.85) +
  scale_fill_viridis_d(begin = 0.8, end = 0) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Gibbs Sampler Iteration",
    y = "Read",
    fill = "Classification"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave(
  z_plot,
  filename = snakemake@output[["z_plot_png"]],
  width = snakemake@params[["plot_width"]],
  height = snakemake@params[["plot_height"]]
)
ggsave(
  z_plot,
  filename = snakemake@output[["z_plot_pdf"]],
  width = snakemake@params[["plot_width"]],
  height = snakemake@params[["plot_height"]]
)