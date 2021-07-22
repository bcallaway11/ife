#' @title ife
#'
#' @description Compute treatment effects in interactive fixed effects models
#'  with a small number of time periods
#'
#' @inheritParams pte:pte
#' @param nife the number of interactive fixed effects to include in the model
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
                anticipation=0,
                alp=0.05,
                biters=100,
                cl=1) {

  res <- pte(yname=yname,
             gname=gname,
             tname=tname,
             idname=idname,
             data=data,
             setup_pte_fun=ife_setup_pte,
             subset_fun=ife_subset,
             attgt_fun=ife_attgt,
             nife=nife,
             xformla=xformla,
             zformla=zformla,
             anticipation=anticipation,
             alp=alp,
             biters=biters,
             cl=cl)

  res
}
