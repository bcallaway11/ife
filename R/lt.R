#' @title linear_trends
#'
#' @description Compute treatment effects in model with linear trends for
#' untreated potential outcomes
#'
#' @inheritParams ife
#'
#' @return \code{ptetools::pte_results} object
#'
#' @export
linear_trends <- function(yname,
                          gname,
                          tname,
                          idname,
                          data,
                          xformla=~1,
                          anticipation=0,
                          cband=TRUE,
                          alp=0.05,
                          boot_type="multiplier",
                          biters=100,
                          cl=1,
                          ...) {

  # set this in order to use `setup_pte` function provided by `pte` package
  required_pre_periods <- 2
  
  res <- pte(yname=yname,
             gname=gname,
             tname=tname,
             idname=idname,
             data=data,
             setup_pte_fun=setup_pte,
             subset_fun=lt_subset,
             attgt_fun=lt_attgt,
             xformla=xformla,
             required_pre_periods=required_pre_periods,
             anticipation=anticipation,
             cband=cband,
             alp=alp,
             boot_type=boot_type,
             biters=biters,
             cl=cl,
             ...)

  res
}
