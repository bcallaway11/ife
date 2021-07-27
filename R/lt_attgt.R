#-----------------------------------------------------------------------------
# Some things to notice:
#
# * We rely on equally spaced periods (though difference doesn't need to be 1)
# * We use all not-yet-treated observations as the comparison group
#    (rather than a fixed never-treated group)
# * This code is a bit of a mess ... 
#-----------------------------------------------------------------------------


#' @title lt_attgt
#' @description Computes estimates of ATT(g,t) using a linear trends model
#'  for untreated potential outcomes
#'
#' @inheritParams ife_attgt
#'
#' @return \code{attgt_if} object
lt_attgt <- function(gt_data,
                      xformla=~1, anticipation=0,  
                      ...) {

  # base period is the first one in this subset of the data
  base.period <- min(gt_data$period)
  tp <- max(gt_data$period)
  this.n <- nrow(gt_data)/3 

  # pick up covariates in the base period (any period ok if
  # they are time constant, otherwise pick up ones in base
  # period
  # for outcome regression, get pre-treatment values
  Xpre <- model.frame(xformla, data=subset(gt_data,period==base.period))
  
  # take difference with respect to base period
  this.data <- gt_data %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(dY_base=(Y-Y[period==base.period])) %>%
    as.data.frame()
  # and drop base period
  this.data <- subset(this.data, period != base.period)

  # account for longer gap between periods for post-treatment periods
  tp_gap <- tp - base.period
  this.data$dY_base <- if_else(this.data$name=="pre", this.data$dY_base*tp_gap, this.data$dY_base)

  pre_data <- subset(this.data, name=="pre")
  post_data <- subset(this.data, name=="post")

  lt_adjusted_Y <- post_data$dY_base - pre_data$dY_base
  post_data$lt_adjusted_Y <- lt_adjusted_Y
  
  # merge outcome and covariate data
  gt_dataX <- cbind.data.frame(post_data[,c("lt_adjusted_Y",
                                            "D",
                                            "id",
                                            "period",
                                            "name")], Xpre)

  # treatment dummy variable
  D <- gt_dataX$D

  # call DRDID functions to make the computations;
  # just like in `did` package
  gt_dataX <- droplevels(gt_dataX)
  attgt <- DRDID::drdid_panel(y1=lt_adjusted_Y,
                              y0=rep(0,length(D)),
                              D=D,
                              covariates=model.matrix(xformla,
                                                      data=gt_dataX),
                              inffunc=TRUE)

  # return attgt
  attgt_if(attgt=attgt$ATT, inf_func=attgt$att.inf.func)
}
