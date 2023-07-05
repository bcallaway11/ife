#' @title staggered_ife
#'
#' @description Compute treatment effects in interactive fixed effects models
#'  with a small number of time periods by exploiting staggered treatment
#'  adoption
#'
#' @inheritParams ife
#' @param nife the number of interactive fixed effects to include in the model
#'
#' @return \code{pte::pte_results} object
#'
#' @export
staggered_ife <- function(yname,
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
  
  # set this in order to use `setup_pte` function provided by `pte` package
  required_pre_periods <- nife+1
  
  res <- pte2(yname=yname,
             gname=gname,
             tname=tname,
             idname=idname,
             data=data,
             setup_pte_fun=staggered_ife_setup_pte,
             subset_fun=ife_subset,
             attgt_fun=staggered_ife_attgt,
             nife=nife,
             xformla=xformla,
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


#' @title staggered_ife2
#'
#' @description Compute treatment effects in interactive fixed effects models
#'  with a small number of time periods by exploiting staggered treatment
#'  adoption.  Unlike `staggered_ife`, this function uses all available 
#'  pre-treatment periods for estimation.  This is the approach taken in 
#'  Callaway and Tsyawo (2023).
#'
#' @inheritParams ife
#' @param nife the number of interactive fixed effects to include in the model
#'
#' @return \code{pte::pte_results} object
#'
#' @export
staggered_ife2 <- function(yname,
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
  
  # set this in order to use `setup_pte` function provided by `pte` package
  required_pre_periods <- nife+1
  
  res <- pte2(yname=yname,
              gname=gname,
              tname=tname,
              idname=idname,
              data=data,
              setup_pte_fun=staggered_ife_setup_data2,
              subset_fun=keep_all_pretreatment_subset,
              attgt_fun=staggered_ife_attgt2,
              nife=nife,
              xformla=xformla,
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
