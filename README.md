
<!-- README.md is generated from README.Rmd. Please edit that file -->

# crrcbcv

<!-- badges: start -->
<!-- badges: end -->

This package offers a user friendly function `crrcbcv` to compute
bias-corrected variances for competing risks regression models using
proportional subdistribution hazards with small-sample clustered data.
Four types of bias correction are included: the MD-type bias correction
by Mancl and DeRouen (2001), the KC-type bias correction by Kauermann
and Carroll (2001), the FG-type bias correction by Fay and Graubard
(2001), and the MBN-type bias correction by Morel, Bokossa, and Neerchal
(2003).

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
#> Loading required package: crrSC
#> Loading required package: survival
#> Loading required package: abind
#> Loading required package: pracma

data(cls)
mod.est = crrc(ftime=cls$T_obs, fstatus=cls$eps, cov1=cls[,c('X_1','X_2')], cluster=cls$I)
crrcbcv(beta=mod.est$coef, ftime=cls$T_obs, fstatus=cls$eps, cov1=cls[,c('X_1','X_2')],
cluster=cls$I, var.type=c('MD','KC','FG','MBN'))
#> $MD
#>            [,1]        [,2]
#> [1,] 0.18712217 0.013304534
#> [2,] 0.01330453 0.007767604
#> 
#> $KC
#>            [,1]        [,2]
#> [1,] 0.14807541 0.155472768
#> [2,] 0.01351378 0.006116421
#> 
#> $FG
#>             [,1]        [,2]
#> [1,] 0.145317028 0.005252785
#> [2,] 0.005252785 0.005509824
#> 
#> $MBN
#>             [,1]        [,2]
#> [1,] 0.133946021 0.004038105
#> [2,] 0.004038105 0.008991263
#> 
#> attr(,"class")
#> [1] "crrcbcv"
```
