SUB_RATES_NEW = [0.005, 0.01, 0.05, 0.1]
PROP_NEW = [0.01, 0.05, 0.1, 0.15, 0.2]
READ_COUNTS = [50, 500, 1000, 2000]

BURNIN = 1000
GIBBS_ITER = 2000

rule all:
    input: 
        expand(
            "plots/pi_dist_plot_{sub_rate_new}.png",
            sub_rate_new = SUB_RATES_NEW
        ),
        "plots/f_dist_plot.png"

rule slample:
    output: 
        f_samples = temp("outputs/f_samples_new_{prop_new}_freads_{read_count}_sub-new_{sub_rate_new}.csv"),
        pi_samples = temp("outputs/pi_g_samples_new_{prop_new}_freads_{read_count}_sub-new_{sub_rate_new}.csv")
    params:
        burnin = BURNIN,
        iterations = GIBBS_ITER
    script: "scripts/slampler.R"

rule plot_pi_g:
    input: 
        pi_samples = expand(
            "outputs/pi_g_samples_new_{prop_new}_freads_{read_count}_sub-new_{{sub_rate_new}}.csv",
            prop_new = PROP_NEW,
            read_count = READ_COUNTS
        )
    output: 
        pi_dist_plot = "plots/pi_dist_plot_{sub_rate_new}.png"
    script: "scripts/plot_pi_g_results.R"

rule plot_f:
    input: 
        pi_samples = expand(
            "outputs/f_samples_new_{prop_new}_freads_{read_count}_sub-new_{sub_rate_new}.csv",
            prop_new = PROP_NEW,
            read_count = READ_COUNTS,
            sub_rate_new = SUB_RATES_NEW
        )
    output: 
        f_dist_plot = "plots/f_dist_plot.png"
    script: "scripts/plot_f_results.R"