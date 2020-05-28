# phyloGLM

<!-- badges: start -->
  [![Travis build status](https://travis-ci.com/ndukler/phyloGLM.svg?branch=master)](https://travis-ci.com/ndukler/phyloGLM)
[![Codecov test coverage](https://codecov.io/gh/ndukler/phyloGLM/branch/master/graph/badge.svg)](https://codecov.io/gh/ndukler/phyloGLM?branch=master)
<!-- badges: end -->

Code to model phylogenetic rate and allelic stationary distribution parameters as function of genomic covariates. Alleles are represented as conditonal probabilities, P(data|allelic state). Allele with only a single value (such as DNA) can be represented by a one-hot encoded vector per-site per-species (e.g. A=[1 0 0 0]). For more detail please see the phyloGLM method
as described and applied in [this](https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msaa073/5807616) paper.

Package can be installed with the following line of code:
`devtools::install_github("ndukler/phyloGLM")`

A prebuilt version of the vignette is published [here](http://rpubs.com/ndukler/https://rpubs.com/ndukler/phyloGLMIntro) but if you want to
rebuild the included vignette to ensure that it is up to date, use:

devtools::install_github("ndukler/phyloGLM", build_vignettes = TRUE)

WARNING: This will take a few minutes.

Acknowledgements:
I would like to thank Nicholas Matzke for useful discussions (https://github.com/nmatzke) and Meabh McCurdy (https://github.com/MeabhMcCurdy) from whom I borrowed and modified code to use expokit routines in c++.
