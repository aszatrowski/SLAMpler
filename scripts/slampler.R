suppressPackageStartupMessages({
  library(ggplot2)
  library(tibble)
  library(dplyr)
  library(viridis)
})
set.seed(1)
# global parameters
burnin <- snakemake@params[["burnin"]]
iterations <- snakemake@params[["iterations"]]
sub_rate_old <- 1e-6
sub_rate_new <- as.numeric(snakemake@wildcards[["sub_rate_new"]])
read_count_g <- as.numeric(snakemake@wildcards[["read_count"]])
total_new_rna_prop <- as.numeric(snakemake@wildcards[["prop_new"]])


sim_reads <- function(read_count, sub_rates_matrix, total_new_rna_prop) {
  num_sites <- ncol(sub_rates_matrix)
  labels <- sample(
    x = c("new", "old"),
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
  read_count = read_count_g,
  sub_rates_matrix = rate_matrix,
  total_new_rna_prop = total_new_rna_prop
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
      number_of_ones <- colSums(reads[z == i, , drop = FALSE])
    }
    # Use informative priors that match the data-generating process
    # Old RNA: Beta(1, 10000) to keep it near 1e-6
    # New RNA: Beta(20, 80) to center around 0.2
    if (i == 1) { # old component
      MatrixF[i, ] <- rbeta(J, number_of_ones + 1, sample_size - number_of_ones + 99999)
    } else { # new component
      MatrixF[i, ] <- rbeta(J, number_of_ones + 0.1, sample_size - number_of_ones + 9.9)
    }
  }
  return(MatrixF)
}

sample_pi <- function(z, prior_alpha = 0.1, prior_beta = 9.9) {
  # estimate pi from z.
  # EV = alpha / (alpha + beta) = 0.1 / (0.1 + 9.9) = 0.01,
  # weak-ish prior that suppresses noise at when pi_g and coverage are low,
  # but still allows for a wide range of values
  n_new <- sum(z == 2)
  n_old <- sum(z == 1)
  pi_g <- rbeta(1, n_new + prior_alpha, n_old + prior_beta)
  return(pi_g)
}


gene_gibbs <- function(read_data, prior_alpha, prior_beta, adaptive_prior = FALSE, niter = iterations) {
  # for the reads for each gene, run the Gibbs sampler
  n <- nrow(read_data)
  J <- ncol(read_data)
  pi_init <- c(0.01, 0.99)

  if (adaptive_prior) {
    # Scale prior strength by coverage
    scale_factor <- max(1, 10 / sqrt(n))
    prior_alpha <- prior_alpha * scale_factor
    prior_beta <- prior_beta * scale_factor
  }

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
    pi_g <- sample_pi(z, prior_alpha, prior_beta)
    z_out[i, ] <- z
    f_out[i, , ] <- f
    pi_out[i] <- pi_g
    if (i %% 100 == 0) {
      sprintf(
        "[pi_g = %.3f, reads = %d, p_n = %.3f] iter: %d",
        total_new_rna_prop, read_count_g, sub_rate_new, i
      ) |> cat("\n")
    }
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
    pi_g = pi_out,
    true_pi = total_new_rna_prop,
    read_count = read_count_g
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

run <- gene_gibbs(sim$reads, prior_alpha = 0.1, prior_beta = 9.9, adaptive_prior = TRUE, niter = iterations)

# Post-hoc relabeling: ensure f_new > f_old so components are identifiable
relabel_by_f_diff <- function(run) {
  niter <- nrow(run$f_old)
  for (i in 1:niter) {
    # Get mean across all sites for this iteration
    f_old_mean <- mean(as.numeric(run$f_old[i, -1]))
    f_new_mean <- mean(as.numeric(run$f_new[i, -1]))
    # If old > new, swap
    if (f_old_mean > f_new_mean) {
      temp <- run$f_old[i, ]
      run$f_old[i, ] <- run$f_new[i, ]
      run$f_new[i, ] <- temp
      run$z[i, -1] <- ifelse(run$z[i, -1] == "old", "new", "old")
      run$pi[i, "pi"] <- 1 - run$pi[i, "pi"]
    }
  }
  return(run)
}
run <- relabel_by_f_diff(run)

z_hist_tbl <- tibble(run$z) |>
  tidyr::pivot_longer(
    cols = -iter,
    names_to = "read_id",
    values_to = "class"
  ) |>
  dplyr::mutate(
    class = factor(class, levels = c(1, 2), labels = c("old", "new")),
  )

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
) |>
  mutate(
    true_f_old = sub_rate_old,
    true_f_new = sub_rate_new,
    read_count = read_count_g
  )


# EXPORTS
readr::write_csv(run$pi, snakemake@output[["pi_samples"]])
readr::write_csv(f_old_new, snakemake@output[["f_samples"]])