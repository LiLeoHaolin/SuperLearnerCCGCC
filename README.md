# Computing Codes for the Paper "Super Learner for Survival Prediction in Case-Cohort and Generalized Case-Cohort Studies"
### Haolin Li, Haibo Zhou, David Couper, and Jianwen Cai

## Description

This repository contains computing codes for the paper "Super Learner for Survival Prediction in Case-Cohort and Generalized Case-Cohort Studies". Please click [here](https://www.google.com) for the full text of the paper. 

## Naming Convention 

### Folders

This repo contains the following folders, corresponding to the simulation codes for reproducing Tables 2 and 3 of the main text. The names and descriptions of the folders are as follows,

* *CC_low_lin* - case-cohort design, low dimensional, linear covariate effect.
* *CC_low_non* - case-cohort design, low dimensional, nonlinear covariate effect.
* *CC_high_lin* - case-cohort design, high dimensional, linear covariate effect.
* *CC_high_non* - case-cohort design, high dimensional, nonlinear covariate effect.
* *GCC_low_lin* - generalized case-cohort design, low dimensional, linear covariate effect.
* *GCC_low_non* - generalized case-cohort design, low dimensional, nonlinear covariate effect.
* *GCC_high_lin* - generalized case-cohort design, high dimensional, linear covariate effect.
* *GCC_high_non* - generalized case-cohort design, high dimensional, nonlinear covariate effect. 

### Files

In each folder, we include the computing codes for all the prediction methods discussed in the paper. The names and descriptions of the files are as follows,

* *dat_gen.r* - the code for data generation.
* *expensive_only.r* - the super learner using the full cohort with inexpensive covariates only.
* *full_cohort.r* - the super learner assuming the full cohort with available expensive covariates for all subjects.
* *naive.r* - the naive super learner, with the same candidate learners as the proposed super learner but without considering design weights in the empirical risk. 
* *proposed.r* - the proposed super learner.
* *SRS.r* - the super learner using a study sample with expensive covariates assembled through an SRS with the same sample size as the case-cohort/generalized case-cohort design.


## References

Li, H., Zhou, H., Couper, D., Cai, J. (2023+). Super Learner for Survival Prediction in Case-Cohort and Generalized Case-Cohort Studies. Manuscript Submitted for Publication.
