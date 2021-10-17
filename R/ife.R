#' @title ife
#'
#' @description Compute treatment effects in interactive fixed effects models
#'  with a small number of time periods
#'
#' @inheritParams ife_attgt
#' @inheritParams pte::pte
#' @param nife the number of interactive fixed effects to include in the model
#' @param cband whether or not to compute a uniform (instead of
#'  pointwise) confidence band
#'
#' @return \code{pte::pte_results} object
#'
#' @export
ife <- function(yname,
                gname,
                tname,
                idname,
                data,
                nife,
                xformla=~1,
                zformla,
                ret_ife_regs=TRUE,
                anticipation=0,
                cband=TRUE,
                alp=0.05,
                boot_type="multiplier",
                biters=100,
                cl=1) {

  # set this in order to use `setup_pte` function provided by `pte` package
  required_pre_periods <- nife+1
  
  res <- pte(yname=yname,
             gname=gname,
             tname=tname,
             idname=idname,
             data=data,
             setup_pte_fun=setup_pte,
             subset_fun=ife_subset,
             attgt_fun=ife_attgt,
             nife=nife,
             xformla=xformla,
             zformla=zformla,
             ret_ife_regs=ret_ife_regs,
             required_pre_periods=required_pre_periods,
             anticipation=anticipation,
             cband=cband,
             alp=alp,
             boot_type=boot_type,
             biters=biters,
             cl=cl)

  res
}
