import numpy as np
import json

SUB_RATES_NEW = np.around(np.arange(0.005, 0.03, 0.005), decimals = 3)
PROP_NEW = np.around(np.arange(0, 1.05, 0.05), decimals=2)
READ_COUNTS = [100, 500, 1000, 2000]

BURNIN = 2500
GIBBS_ITER = 5000
if BURNIN >= GIBBS_ITER * 0.6:
    raise ValueError("Too many iterations discarded as burn-in")

PLOT_WIDTH = 8
PLOT_HEIGHT = 4

wildcards_to_plot = dict(
    prop_new = '0.1',
    read_count = '500',
    sub_rate_new = '0.01'
)

localrules: plot_pi_g, plot_f, expected_subs, plot_z

rule all:
    input: 
        expand(
            "plots/pi_dist_plot_{sub_rate_new}.png",
            sub_rate_new = SUB_RATES_NEW
        ),
        expand(
            "plots/pi_dist_plot_{sub_rate_new}.pdf",
            sub_rate_new = SUB_RATES_NEW
        ),
        "plots/f_dist_plot.pdf",
        "plots/f_dist_plot.png",
        "plots/expected_subs.pdf",
        "plots/expected_subs.png",
        "plots/z_classifications.pdf",
        "plots/z_classifications.png"

rule slample:
    output: 
        f_samples = "samples/f_samples_new_{prop_new}_reads_{read_count}_sub-new_{sub_rate_new}.csv",
        pi_samples = "samples/pi_g_samples_new_{prop_new}_reads_{read_count}_sub-new_{sub_rate_new}.csv",
        z_samples = "samples/z_samples_new_{prop_new}_reads_{read_count}_sub-new_{sub_rate_new}.csv"
    params:
        burnin = BURNIN,
        iterations = GIBBS_ITER,
        prior_alpha = 0.8,
        prior_beta = 0.9,
        wildcards_to_plot = json.dumps(wildcards_to_plot)
    conda: 'env.yaml'
    resources:
        runtime = 25,
        mem = "4GB"
    benchmark:
        "benchmarks/bench_{prop_new}_reads_{read_count}_sub-new_{sub_rate_new}.txt"
    script: "scripts/slampler.R"

rule plot_pi_g:
    input: 
        pi_samples = expand(
            "samples/pi_g_samples_new_{prop_new}_reads_{read_count}_sub-new_{{sub_rate_new}}.csv",
            prop_new = PROP_NEW,
            read_count = READ_COUNTS
        )
    output: 
        pi_dist_plot_pdf = "plots/pi_dist_plot_{sub_rate_new}.pdf",
        pi_dist_plot_png = "plots/pi_dist_plot_{sub_rate_new}.png"
    params:
        plot_height = PLOT_HEIGHT,
        plot_width = PLOT_WIDTH
    conda: 'env.yaml'
    resources:
        runtime = 5,
        mem = "2GB"
    script: "scripts/plot_pi_g_results.R"

rule plot_f:
    input: 
        pi_samples = expand(
            "samples/f_samples_new_{prop_new}_reads_{read_count}_sub-new_{sub_rate_new}.csv",
            prop_new = PROP_NEW,
            read_count = READ_COUNTS,
            sub_rate_new = SUB_RATES_NEW
        )
    output: 
        f_dist_plot_pdf = "plots/f_dist_plot.pdf",
        f_dist_plot_png = "plots/f_dist_plot.png"
    params:
        plot_height = PLOT_HEIGHT,
        plot_width = PLOT_WIDTH
    conda: 'env.yaml'
    resources:
        runtime = 5,
        mem = "2GB"
    script: "scripts/plot_f_results.R"

rule expected_subs:
    """
    Plot distribution of number of substitutions on a 150bp read given a substitution rate.
    """
    output:
        binom_plot_pdf = "plots/expected_subs.pdf",
        binom_plot_png = "plots/expected_subs.png"
    params: 
        sub_rates = SUB_RATES_NEW,
        plot_height = PLOT_HEIGHT,
        plot_width = PLOT_WIDTH
    conda: 'env.yaml'
    script: "scripts/expected_subs.R"

rule plot_z:
    input: "samples/z_samples_new_0.1_reads_500_sub-new_0.01.csv"
    output: 
        z_plot_pdf = "plots/z_classifications.pdf",
        z_plot_png = "plots/z_classifications.png"
    params:
        plot_height = PLOT_HEIGHT,
        plot_width = PLOT_WIDTH
    conda: 'env.yaml'
    resources:
        runtime = 5,
        mem = "2GB"
    script: "scripts/plot_z_results.R"