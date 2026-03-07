SUB_RATES_NEW = [0.005, 0.010, 0.015, 0.020]
PROP_NEW = [0.01, 0.1, 0.2]
READ_COUNTS = [50, 500, 1000, 2000]

BURNIN = 500
GIBBS_ITER = 1000
if BURNIN >= GIBBS_ITER * 0.6:
    raise ValueError("Too many iterations discarded as burn-in")

localrules: plot_pi_g, plot_f

rule all:
    input: 
        expand(
            "plots/pi_dist_plot_{sub_rate_new}.png",
            sub_rate_new = SUB_RATES_NEW
        ),
        "plots/f_dist_plot.png"

rule slample:
    output: 
        f_samples = "samples/f_samples_new_{prop_new}_reads_{read_count}_sub-new_{sub_rate_new}.csv",
        pi_samples = "samples/pi_g_samples_new_{prop_new}_reads_{read_count}_sub-new_{sub_rate_new}.csv"
    params:
        burnin = BURNIN,
        iterations = GIBBS_ITER
    conda: 'env.yaml'
    resources:
        runtime = 15,
        mem = "4GB"
    script: "scripts/slampler.R"

rule plot_pi_g:
    input: 
        pi_samples = expand(
            "samples/pi_g_samples_new_{prop_new}_reads_{read_count}_sub-new_{{sub_rate_new}}.csv",
            prop_new = PROP_NEW,
            read_count = READ_COUNTS
        )
    output: 
        pi_dist_plot = "plots/pi_dist_plot_{sub_rate_new}.png"
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
        f_dist_plot = "plots/f_dist_plot.png"
    conda: 'env.yaml'
    resources:
        runtime = 5,
        mem = "2GB"
    script: "scripts/plot_f_results.R"