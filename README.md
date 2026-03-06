# HGEN 486 Final Project: The Gibbs SLAMpler
*Austin Szatrowski*

A Gibbs sampling approach for analysis of SLAM-seq data.

## Ideas:
- [x] (1) snakemake wrapper
- [ ] (2) Performance axes
    - [x] (2.1) over a range of $\pi_g$
        - [x] plots of distributions vs true—aggregation mechanism via snakemake and then plot performance
        - [ ] get the prior calibration right
    - [ ] (2.2) read counts
        - [x] run
        - [x] calibration curves
        - [ ] strengthen prior for low read counts—little bit of empirical Bayes 
    - [ ] (2.3) substitution rates (old and new, e.g. sequencing error and 4sU incorporation) (experimental quality)
        - [ ] accuracy heatmap f(x = old, y = new)
        - [ ] or just separate sets of calibration curves
- [ ] move to cluster
- [ ] biologically informed priors
    - if (likely activated, via gene annotation):
        - **a module to interface with GSEA**
        - prior on new RNA ++
- [ ] MCMC diagnostics with `coda`
- [ ] Experiment with priors
- [ ] Derivations on poster (maybe)
- [ ] Build a BAM file interface, in particular one that is sensitive to alignment quality
- [ ] identify reads that really don't fit the pattern?? does that make statistical sense??
- [ ] study substitution rates both by gene and by site
- [ ] think more about optimizing $\pi_g$ more than $z$

## Notes for poster:
- this is a sparse sampling problem where informative priors will be most useful
`rcartocolor::Geyser`
`nationalparkcolors::MtRainier`
` nationalparkcolors::GrandTeton`
`ltc::alger`
`ltc::franscoise`