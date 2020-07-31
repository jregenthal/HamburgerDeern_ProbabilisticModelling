# Project for Probabilistic Modelling by Hamburger Deern

This project evaluates several methods to compare Bayesian models for claim count data. It is based on the paper "Bayesian Assessment of the Distribution of Insurance Claim Counts Using Reversible Jump MCMC" by Ntzoufras et al. (2005). 

The files are structured in six parts accourding to their prefixes:

* 00_: Files regarding the original data and data simulation
* 01_: R codes for the implementation of the initial models in JAGS
* 02_: R codes for the RJMCMC from scratch
* 03_: R codes for the post-processing approach
* 04_: R codes for the bridge sampling approach including Bayes' factor and DIC
* 05_: R codes for the WAIC