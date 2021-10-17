#' @title cov_pre_test
#'
#' @description Pre-test for time invariant effects of covariates in an
#'  interactive fixed effects model for untreated potential outcomes
#'
#' @inheritParams ife
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
                         anticipation=0) {


  required_pre_periods <- nife+1
  
  ptep <- setup_cov_pre_test(yname=yname,
                             gname=gname,
                             tname=tname,
                             idname=idname,
                             data=data,
                             anticipation=anticipation,
                             cband=FALSE,
                             alp=.05,
                             boot_type="multiplier",
                             biters=100,
                             cl=1)
  
  res <- compute.pte(ptep=ptep,
                     subset_fun=cov_pre_test_subset,
                     attgt_fun=cov_pre_test_attgt,
                     nife=nife,
                     xformla=xformla,
                     anticipation=anticipation,
                     ret_ife_regs=TRUE)
  
  out_regs <- res$extra_gt_returns
  keepers <- sapply(out_regs, function(or) {
    (!is.na(or$extra_gt_returns)) & (or$group - or$time.period > anticipation) 
  })


  # groups and time periods back to originals
  out_regs <- out_regs[keepers]
  original_t <- unique(data[,tname])
  out_regs <- lapply(out_regs, function(or) {
    or$group <- t2orig(or$group, original_t)
    or$time.period <- t2orig(or$time.period, original_t)
    or
  })

  class(out_regs) <- "cov_pretest_results"
  out_regs
}


#' @title summary.cov_pretest_results
#'
#' @param object a cov_pretest_results object
#' @param ... other arguments
#'
#' @export
summary.cov_pretest_results <- function(object, ...) {
  gt_ests <- lapply(object, function(ob) ob$extra_gt_returns$ife_reg)
  names(gt_ests) <- sapply(object, function(ob) {
    paste0("g:",ob$group, ",t:", ob$time.period)
  })
  modelsummary(gt_ests, stars=c("*"=.05), output="markdown")
}
