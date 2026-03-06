# SUB_RATES_NEW = [1e-4, 0.01, 0.05, 0.1, 0.15, 0.2]
# PROP_NEW = [0, 1e-3, 0.01, 0.05, 0.1, 0.15, 0.2]
PROP_NEW = [0.01, 0.05, 0.1, 0.15, 0.2]
READ_COUNTS = [50, 500, 1000, 2000]

BURNIN = 1000
GIBBS_ITER = 2000

rule all:
    input: 
        expand(
            "outputs/f_samples_new_{prop_new}_reads_{read_count}.csv",
            prop_new = PROP_NEW,
            read_count = READ_COUNTS,
        ),
        expand(
            "outputs/pi_g_samples_new_{prop_new}_reads_{read_count}.csv",
            prop_new = PROP_NEW,
            read_count = READ_COUNTS,
        ),
        "plots/pi_dist_plot.png"

rule slample:
    output: 
        f_samples = "outputs/f_samples_new_{prop_new}_reads_{read_count}.csv",
        pi_samples = "outputs/pi_g_samples_new_{prop_new}_reads_{read_count}.csv",
    params:
        burnin = BURNIN,
        iterations = GIBBS_ITER
    script: "scripts/slampler.R"

rule plot_pi_g:
    input: 
        pi_samples = expand(
            "outputs/pi_g_samples_new_{prop_new}_reads_{read_count}.csv",
            prop_new = PROP_NEW,
            read_count = READ_COUNTS
        ),
    output: 
        pi_dist_plot = "plots/pi_dist_plot.png"
    script: "scripts/plot_pi_g_results.R"