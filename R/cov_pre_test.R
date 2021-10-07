#' @title cov_pre_test
#'
#' @description Pre-test for time invariant effects of covariates in an
#'  interactive fixed effects model for untreated potential outcomes
#'
#' @inheritParams pte:pte
#' @param nife the number of interactive fixed effects to include in the model
#'
#' @return \code{pte::pte_results} object
#'
#' @export
cov_pre_test <- function(yname,
                         gname,
                         tname,
                         idname,
                         data,
                         nife,
                         xformla=~1,
                         ret_ife_regs=TRUE,
                         anticipation=0,
                         cband=TRUE,
                         alp=0.05,
                         boot_type="multiplier",
                         biters=100,
                         cl=1) {


  required_pre_periods <- nife+1
  
  ptep <- setup_cov_pre_test(yname=yname,
                             gname=gname,
                             tname=tname,
                             idname=idname,
                             data=data,
                             cband=cband,
                             alp=alp,
                             boot_type=boot_type,
                             biters=biters,
                             cl=cl)
  
  res <- compute.pte(ptep=ptep,
                     subset_fun=cov_pre_test_subset,
                     attgt_fun=cov_pre_test_attgt,
                     nife=nife,
                     xformla=xformla,
                     anticipation=anticipation,
                     ret_ife_regs=TRUE)
  
  out_regs <- res$extra_gt_returns
  keepers <- sapply(out_regs, function(or) !is.na(or$extra_gt_returns))

  out_regs[keepers]
}
