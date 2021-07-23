
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<!-- badges: end -->

# Interactive Fixed Effects (ife) Package

The `ife` package contains code to estimate treatment effects in a setup
where a researcher has access to panel data (or, hopefully in the near
future, repeated cross sections data) and where untreated potential
outcomes are generated by an interactive fixed effects model.

The package is not currently available on CRAN, but it can be the
development version of the package can be installed from github by

``` r
# install.packages("devtools")
devtools::install_github("bcallaway11/ife")
```

## Example

Next, we provide a brief example using the application from Callaway and
Karami (2021).

``` r
res <- ife(yname="earn",
           gname="first.displaced",
           tname="year",
           idname="id",
           data=job_displacement_data,
           nife=1,
           xformla=~EDUC + race + gender,
           zformla=~EDUC + race + gender + afqt,
           anticipation=1,
           alp=0.1,
           biters=1000)

summary(res)
#> 
#> Overall ATT:  
#>        ATT    Std. Error     [ 90%  Conf. Int.]  
#>  -3913.851      1192.724  -5875.707   -1951.994 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time   Estimate Std. Error     [90%  Conf. Band]  
#>          -6   -63.8360   939.6738 -2209.276   2081.6039  
#>          -4   -63.4326   699.6486 -1660.853   1533.9877  
#>          -2  -664.5203   617.0093 -2073.261    744.2200  
#>           0 -3848.8739  1051.8789 -6250.498  -1447.2502 *
#>           2 -4179.1319  1492.4485 -7586.653   -771.6106 *
#>           4 -1454.0791  2640.8284 -7483.553   4575.3947  
#> ---
#> Signif. codes: `*' confidence band does not cover 0
ggpte(res)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />
