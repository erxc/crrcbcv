
<!-- README.md is generated from README.Rmd. Please edit that file -->

# crrcbcv

<!-- badges: start -->
<!-- badges: end -->

This package offers a user friendly function 'crrcbcv' to compute bias-corrected variances for competing risks regression models using proportional subdistribution hazards with small-sample clustered data. Four types of bias correction are included: the MD-type bias correction by Mancl and DeRouen (2001), the KC-type bias correction by Kauermann and Carroll (2001), the FG-type bias correction by Fay and Graubard (2001), and the MBN-type bias correction by Morel, Bokossa, and Neerchal (2003).

## Installation

You can install the released version of crrcbcv from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("crrcbcv")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(crrcbcv)
#> Loading required package: survival
#> Warning: package 'survival' was built under R version 4.1.1
#> Loading required package: crrSC
#> Loading required package: abind
#> Loading required package: pracma
## basic example code

data(cls)
mod.est = crrc(ftime=cls$T_obs, fstatus=cls$eps, cov1=cls[,c('X_1','X_2')], cluster=cls$I)
crrcbcv(beta=mod.est$coef, ftime=cls$T_obs, fstatus=cls$eps, cov1=cls[,c('X_1','X_2')],
cluster=cls$I, var.type=c('MD','KC','FG','MBN'))
```
