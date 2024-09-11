
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
#>  -3976.881      1190.287  -5934.728   -2019.033 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time   Estimate Std. Error     [90%  Conf. Band]  
#>          -6   -15.0413   811.6706 -1350.121   1320.0380  
#>          -4  -131.3128   703.7811 -1288.930   1026.3041  
#>          -2  -664.5203   580.6654 -1619.630    290.5892  
#>           0 -3946.5861  1036.9978 -5652.296  -2240.8765 *
#>           2 -4281.1615  1355.4629 -6510.699  -2051.6234 *
#>           4 -1454.0791  2545.9788 -5641.842   2733.6833  
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
#>  -3463.467      1324.384  -5641.884   -1285.049 *
#> 
#> 
#> Dynamic Effects:
#>  Event Time   Estimate Std. Error     [90%  Conf. Band]  
#>          -6   -17.6954   911.3062 -1516.661   1481.2698  
#>          -4   -69.4965   765.2249 -1328.179   1189.1864  
#>          -2  -586.0360   677.8071 -1700.929    528.8575  
#>           0 -3599.9237  1125.4559 -5451.134  -1748.7135 *
#>           2 -3628.0465  1716.6073 -6451.614   -804.4788 *
#>           4   326.6318  3205.9479 -4946.683   5599.9469  
#> ---
#> Signif. codes: `*' confidence band does not cover 0
ggpte(lt_res) + ylim(c(-7000,7000))
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />
