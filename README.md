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
        - [x] adjust prior on $f$ downward
        - [ ] weaken prior for low read counts
        - [ ] eBayes: get $\pi_g$ MLE from the data
    - [ ] (2.3) substitution rates (old and new, e.g. sequencing error and 4sU incorporation) (experimental quality)
        - [ ] accuracy heatmap f(x = old, y = new)
        - [ ] or just separate sets of calibration curves
- [x] move to cluster
- [ ] biologically informed priors
    - if (likely activated, via gene annotation):
        - **a module to interface with GSEA**
        - prior on new RNA ++
- [ ] Build a BAM file interface, in particular one that is sensitive to alignment quality
- [ ] identify reads that really don't fit the pattern?? does that make statistical sense??
- [ ] MCMC diagnostics with `coda`
- [ ] Derivations on poster (maybe)

## Notes for poster:
- this is a sparse sampling problem where informative priors will be most useful
`rcartocolor::Geyser`
`nationalparkcolors::MtRainier`
` nationalparkcolors::GrandTeton`
`ltc::alger`
`ltc::franscoise`