SUB_RATES_NEW = [1e-4, 0.01, 0.05, 0.1, 0.15, 0.2]
PROP_NEW = [0, 1e-3, 0.01, 0.05, 0.1, 0.15, 0.2]

BURNIN = 1000
GIBBS_ITER = 2000
rule all:
    input: 
        expand(
            "outputs/f_samples_new_{prop_new}.csv",
            prop_new = PROP_NEW
        ),
        expand(
            "outputs/pi_g_samples_new_{prop_new}.csv",
            prop_new = PROP_NEW
        ),
        expand(
            "plots/{plot}_plot_new_{prop_new}.png",
            prop_new = PROP_NEW,
            plot = ['pi']
        ),
        "plots/pi_dist_plot.png"

rule slample:
    output: 
        f_samples = "outputs/f_samples_new_{prop_new}.csv",
        pi_samples = "outputs/pi_g_samples_new_{prop_new}.csv",
        pi_plot = "plots/pi_plot_new_{prop_new}.png",
    params:
        burnin = BURNIN,
        iterations = GIBBS_ITER
    script: "scripts/slampler.R"

rule plot_pi_g:
    input: 
        pi_samples = expand(
            "outputs/pi_g_samples_new_{prop_new}.csv",
            prop_new = PROP_NEW
        ),
    output: 
        pi_dist_plot = "plots/pi_dist_plot.png"
    script: "scripts/plot_pi_g_results.R"