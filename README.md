
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

<!-- badges: end -->

# Interactive Fixed Effects (ife) Package <img src="man/figures/logo.png" align="right" height="139" alt="" />

The `ife` package contains code to estimate treatment effects in a setup
where a researcher has access to panel data (or, hopefully in the near
future, repeated cross sections data) and where untreated potential
outcomes are generated by an interactive fixed effects model.

The package is not currently available on CRAN, but the development
version of the package can be installed from github by

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
           ret_ife_regs=TRUE,
           anticipation=1,
           cband=FALSE,
           alp=0.10,
           boot_type="multiplier",
           biters=1000,
           cl=10)

summary(res)
#> 
#> Overall ATT:  
#>        ATT    Std. Error     [ 90%  Conf. Int.]  
#>  -3976.881      1141.811  -5854.993   -2098.769 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time   Estimate Std. Error     [90%  Conf. Band]  
#>          -6   -15.0413   791.5992 -1317.106   1287.0234  
#>          -4  -131.3128   701.6462 -1285.418   1022.7925  
#>          -2  -664.5203   580.2852 -1619.005    289.9638  
#>           0 -3946.5861   986.6150 -5569.423  -2323.7488 *
#>           2 -4281.1615  1537.6912 -6810.438  -1751.8846 *
#>           4 -1454.0791  2591.7845 -5717.185   2809.0271  
#> ---
#> Signif. codes: `*' confidence band does not cover 0
ggpte(res) + ylim(c(-7000,7000))
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

We also have some code for running *individual-specific linear trends
models*. These are a special case of the interactive fixed effects
models that we consider in the paper, but where the factors
![F\_t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_t
"F_t") are restricted to be equal to
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t
"t"). We mostly argue *against* these sorts of models in the paper, but
one advantage is that they do not require any restrictions/assumptions
about finding a covariate whose effects do not change over time.

This code also implements a version of linear trends that is specific to
untreated potential outcomes. Presumably, many of the same criticisms
(and perhaps more actually) in recent papers about implementing DID with
a two-way fixed effects regression likely apply when one includes
individual-specific linear trends in the same sort of specification. The
code we provide here circumvents those issues.

``` r
lt_res <- linear_trends(yname="earn",
                        gname="first.displaced",
                        tname="year",
                        idname="id",
                        data=job_displacement_data,
                        xformla=~EDUC + race + gender,
                        anticipation=1,
                        cband=FALSE,
                        alp=0.10,
                        boot_type="multiplier",
                        biters=1000,
                        cl=10)

summary(lt_res)
#> 
#> Overall ATT:  
#>        ATT    Std. Error     [ 90%  Conf. Int.]  
#>  -3463.467      1350.412  -5684.697   -1242.236 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time   Estimate Std. Error     [90%  Conf. Band]  
#>          -6   -17.6954   823.6747 -1372.520   1337.1290  
#>          -4   -69.4965   702.3175 -1224.706   1085.7130  
#>          -2  -586.0360   654.9993 -1663.414    491.3420  
#>           0 -3599.9237  1205.0032 -5581.978  -1617.8697 *
#>           2 -3628.0465  1808.2979 -6602.432   -653.6612 *
#>           4   326.6318  3178.3581 -4901.302   5554.5656  
#> ---
#> Signif. codes: `*' confidence band does not cover 0
ggpte(lt_res) + ylim(c(-7000,7000))
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />
