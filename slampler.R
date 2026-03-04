library(ggplot2)
library(tibble)
library(dplyr)
library(viridis)
set.seed(1)

# global parameters
burnin <- 500
iterations <- 1000
sub_rate_old <- 1e-6
sub_rate_new <- 0.2


sim_reads <- function(read_count, sub_rates_matrix, total_new_rna_prop) {
  num_sites <- ncol(sub_rates_matrix)
  labels <- sample(
    c("old", "new"),
    prob = c(total_new_rna_prop, 1 - total_new_rna_prop),
    size = read_count,
    replace = TRUE
  )
  reads <- matrix(0, read_count, num_sites)
  for (i in seq_len(read_count)) {
    reads[i, ] <- rbinom(num_sites, rep(1, num_sites), sub_rates_matrix[labels[i], ])
  }
  return(list(reads = reads, labels = labels))
}

sub_rate_matrix <- function(sub_rate_old, sub_rate_new, read_length) {
  # given T>C substitution rates in old and new RNA, build a 2 x 150bp matrix of per-site rates
  sub_rates <- matrix(
    c(
      rep(sub_rate_old, read_length),
      rep(sub_rate_new, read_length)
    ),
    nrow = 2,
    byrow = TRUE
  )
  rownames(sub_rates) <- c("old", "new")
  return(sub_rates)
}
rate_matrix <- sub_rate_matrix(
  sub_rate_old = sub_rate_old,
  sub_rate_new = sub_rate_new,
  read_length = 150
)
sim <- sim_reads(
  read_count = 1000,
  sub_rates_matrix = rate_matrix,
  total_new_rna_prop = 0.5
)

normalize <- function(x) {
  return(x / sum(x))
}
loglik_read_rates <- function(read, f) {
  # read: vector of length J (sites), binary 0/1
  # f: 2 x J matrix (row 1 = old rates, row 2 = new rates)
  # returns: vector of length 2 (loglik of this read under old, under new)
  # for each component (old/new), calculate log P(read | rates)
  # loglik = sum over sites of: x_j * log(f) + (1-x_j) * log(1-f)
  loglik <- rep(NA, 2)
  for (k in 1:2) {
    loglik[k] <- sum(read * log(f[k, ]) + (1 - read) * log(1 - f[k, ]))
  }
  return(loglik)
}

sample_z <- function(reads, f, pi_g) {
  # reads: N x J matrix of reads (rows = reads, cols = sites)
  # f: 2 x J matrix of conversion rates (row 1 = old, row 2 = new)
  # pi_g: scalar, current proportion of new RNA
  # returns: vector of length N with assignments (1 = old, 2 = new)
  N <- nrow(reads)
  # calculate log-likelihood of each read under each component
  loglik_matrix <- apply(
    X = reads,
    MARGIN = 1,  # iterate over reads
    FUN = loglik_read_rates,
    f = f
  )
  # loglik_matrix is now 2 x N (row 1 = loglik under old, row 2 = loglik under new)
  # add log prior (mixing proportions)
  log_pi <- log(c(1 - pi_g, pi_g))
  loglik_matrix <- loglik_matrix + log_pi  # broadcasting across columns
  # convert to probabilities (exponentiate and normalize)
  lik_matrix <- exp(loglik_matrix)
  p_z_given_read <- apply(
    X = lik_matrix,
    MARGIN = 2,  # normalize each column (each read)
    FUN = normalize
  )
  # p_z_given_read is now 2 x N, with columns summing to 1
  # sample assignments
  z <- numeric(N)
  for (i in 1:N) {
    z[i] <- sample(
      x = 1:2,
      size = 1,
      prob = p_z_given_read[, i]
    )
  }
  return(z)
}

sample_f <- function(reads, z) {
  J <- ncol(reads) # number of sites
  MatrixF <- matrix(NA, nrow = 2, ncol = J) # allele freq matrix, 2 x J
  for (i in seq_len(2)) {
    sample_size <- sum(z == i)
    if (sample_size == 0) {
      number_of_ones <- rep(0, J)
    } else {
      # when only one read is assigned to a component the subsetting
      # returns a vector, and colSums() will complain that the input
      # has fewer than two dimensions.  force a matrix with drop = FALSE.
      number_of_ones <- colSums(reads[z == i, , drop = FALSE])
    }
    MatrixF[i, ] <- rbeta(J, number_of_ones + 1, sample_size - number_of_ones + 1)
  }
  return(MatrixF)
}

sample_pi <- function(z, prior_alpha = 1, prior_beta = 5) {
  # estimate pi from z.
  # if prior_beta >> prior_alpha, put a stronger prior on new RNA being rare
  n_new <- sum(z == 2)
  n_old <- sum(z == 1)
  pi_g <- rbeta(1, n_new + prior_alpha, n_old + prior_beta)
  return(pi_g)
}

gene_gibbs <- function(read_data, niter = iterations) {
  # for the reads for each gene, run the Gibbs sampler
  n <- nrow(read_data)
  J <- ncol(read_data)
  pi_init <- c(0.1, 0.9)

  z_out <- matrix(NA, nrow = niter, ncol = n)
  pi_out <- rep(NA, niter)
  f_out <- array(0, dim = c(niter, 2, J))

  z <- sample(x = 2, size = n, replace = TRUE, prob = pi_init)
  z_out[1, ] <- z

  f <- sample_f(read_data, z)
  f_out[1, , ] <- f

  pi_g <- sample_pi(z)
  pi_out[1] <- pi_g

  for (i in 2:niter) {
    f <- sample_f(read_data, z)
    z <- sample_z(read_data, f, pi_g)
    pi_g <- sample_pi(z, prior_alpha = 1, prior_beta = 5)
    z_out[i, ] <- z
    f_out[i, , ] <- f
    pi_out[i] <- pi_g
    print(paste("Iteration", i, "complete"))
  }
  z_out_df <- data.frame(
    iter = seq_len(nrow(z_out)),
    z = z_out
  )[burnin:nrow(z_out), ]

  f_old_out <- f_out[, 1, ]  # old RNA rates
  f_new_out <- f_out[, 2, ]  # new RNA rates

  f_old_out_df <- data.frame(
    iter = seq_len(nrow(f_old_out)),
    f = f_old_out
  )[burnin:nrow(f_old_out), ]
  f_new_out_df <- data.frame(
    iter = seq_len(nrow(f_new_out)),
    f = f_new_out
  )[burnin:nrow(f_new_out), ]

  pi_out_df <- data.frame(
    iter = seq_along(pi_out),
    pi = pi_out
  )[burnin:length(pi_out), ]

  return(
    list(
      z = z_out_df,
      f_old = f_old_out_df,
      f_new = f_new_out_df,
      pi = pi_out_df
    )
  )
}

run <- gene_gibbs(sim$reads)

z_hist_tbl <- tibble(run$z) |>
  tidyr::pivot_longer(
    cols = -iter,
    names_to = "read_id",
    values_to = "class"
  ) |>
  dplyr::mutate(
    class = factor(class, levels = c(1, 2), labels = c("old", "new")),
  )

classification_plot <- ggplot(z_hist_tbl, aes(x = iter, y = read_id)) +
  geom_raster(aes(fill = class)) +
  scale_fill_viridis_d(begin = 0.25, end = 0.75) +
  labs(
    x = "Gibbs Sampler Iteration",
    y = "Read",
    fill = "Population Estimate"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  )
ggsave(
  classification_plot,
  filename = "plots/classification_plot.png",
  width = 6,
  height = 8
)

pi_g_plot <- ggplot(run$pi, aes(x = iter, y = pi)) +
  geom_line(
    color = viridis(1, alpha = 1, begin = 0.5, end = 1, direction = 1, option = "D")
  ) +
  labs(
    x = "Gibbs Sampler Iteration",
    y = "Estimated Proportion of New RNA"
  ) +
  theme_bw()

ggsave(
  pi_g_plot,
  filename = "plots/pi_g_plot.png",
  width = 6,
  height = 8
)

pi_g_dist_plot <- ggplot(run$pi, aes(x = pi)) +
  geom_histogram(
    bins = 30,
    fill = viridis(1, alpha = 1, begin = 0.5, end = 1, direction = 1, option = "D"),
  ) +
  labs(
    x = "Estimated Proportion of New RNA",
    y = "Frequency"
  ) +
  theme_bw()
ggsave(
  pi_g_dist_plot,
  filename = "plots/pi_g_dist_plot.png",
  width = 6,
  height = 8
)

readr::write_csv(run$pi, "outputs/pi_g_samples.csv")

pivot_f <- function(f_df) {
  f_long <- f_df |>
    tidyr::pivot_longer(
      cols = -iter,
      names_to = "site",
      values_to = "f"
    ) |>
    group_by(iter) |>
    summarise(f = mean(f))
  return(f_long)
}

f_old_new <- inner_join(
  pivot_f(run$f_old) |> rename(f_old = f),
  pivot_f(run$f_new) |> rename(f_new = f),
  by = "iter"
)
f_old_new_long <- f_old_new |>
  tidyr::pivot_longer(
    cols = c(f_old, f_new),
    names_to = "population",
    values_to = "f"
  ) |>
  mutate(
    population = factor(population, levels = c("f_old", "f_new"), labels = c("old", "new"))
  )

f_plot <- ggplot(f_old_new_long, aes(x = iter, y = f)) +
  geom_line(aes(color = population)) +
  scale_color_viridis_d(begin = 0.25, end = 0.75) +
  geom_hline(yintercept = c(sub_rate_new, sub_rate_old), linetype = "dashed", color = "gray") +
  labs(
    x = "Gibbs Sampler Iteration",
    y = "Estimated Substitution Rate in new RNA"
  ) +
  theme_bw()
ggsave(
  f_plot,
  filename = "plots/f_plot.png",
  width = 6,
  height = 8
)

readr::write_csv(f_old_new, "outputs/f_samples.csv")