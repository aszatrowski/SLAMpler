# HGEN 486 Final Project: The Gibbs SLAMpler
*Austin Szatrowski*

A Gibbs sampling approach for analysis of SLAM-seq data.

See the [poster](26-03-09_Szatrowski_FinalPoster.pdf) for details.

## Todo:
- [x] (1) snakemake wrapper
- [ ] (2) Performance axes
    - [x] (2.1) over a range of $\pi_g$
        - [x] plots of distributions vs true—aggregation mechanism via snakemake and then plot performance
        - [ ] get the prior calibration right
    - [x] (2.2) read counts
        - [x] run
        - [x] calibration curves
        - [x] adjust prior on $f$ downward
        - [x] weaken prior for low read counts—now fixed.
            - [x] bump up iterations
            - [x] add a few more values
- [x] move to cluster
- [x] Set prior to $\text{Beta}(0.8, 0.9)$ to match empirical distribution
- [x] POSTER
    - [x] initial set of figures
    - [x] stability of classification
- [ ] eBayes: get $\pi_g$ MLE from the data
- [ ] MCMC diagnostics with `coda`
    - [ ] genes with more reads will stabilize faster, so can use fewer iterations to achieve a given degree of precision
- [ ] implement the GRAND-SLAM algorithm in R and compare the results
- [ ] Build a BAM file interface, in particular one that is sensitive to alignment quality and propagate that uncertainty
- [ ] identify reads that really don't fit the pattern?? does that make statistical sense??