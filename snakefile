SUB_RATES_NEW = [0.01, 0.05, 0.1, 0.2]
PROP_NEW = [0.01, 0.05, 0.1, 0.2]
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
            plot = ['classification', 'pi_dist', 'pi']
        )

rule slample:
    output: 
        f_samples = "outputs/f_samples_new_{prop_new}.csv",
        pi_samples = "outputs/pi_g_samples_new_{prop_new}.csv",
        pi_plot = "plots/pi_plot_new_{prop_new}.png",
        classification_plot = "plots/classification_plot_new_{prop_new}.png",
        pi_dist_plot = "plots/pi_dist_plot_new_{prop_new}.png"
    script: "scripts/slampler.R"