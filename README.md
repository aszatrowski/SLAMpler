# HGEN 486 Final Project: The Gibbs SLAMpler
*Austin Szatrowski*

A Gibbs sampling approach for analysis of SLAM-seq data.

## Ideas:
- [ ] (1) snakemake wrapper
- [ ] (2) evaluate performance over a range of true $\pi_g$ and substitution rates $p_n$ and $p_o$
    - [ ] $F_1$ score or similar for binary classification accuracy
    - [ ] plots of distributions vs true—aggregation mechanism via snakemake and then plot performance
    - [ ] also over a range of substitution rates; hold sequencing errors constant
- [ ] (3) allow user to specify priors for incorporation rates
- [x] good output for estimating substitution rates, since that's really the innovation we're trying to achieve here
    - both on same plot
- [ ] MCMC diagnostics with `coda`
- [ ] gene annotation priors (as discussed w JN)
- [ ] Experiment with priors
- [ ] Derivations on poster (maybe)
- [ ] Build a BAM file interface, in particular one that is sensitive to alignment quality
- [ ] identify reads that really don't fit the pattern?? does that make statistical sense??
- [ ] study substitution rates both by gene and by site
- [ ] think more about optimizing $\pi_g$ more than $z$

`rcartocolor::Geyser`
`nationalparkcolors::MtRainier`
` nationalparkcolors::GrandTeton`
`ltc::alger`
`ltc::franscoise`