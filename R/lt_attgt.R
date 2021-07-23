#-----------------------------------------------------------------------------
# Some things to notice:
#
# * We rely on equally spaced periods (though difference doesn't need to be 1)
# * We use all not-yet-treated observations as the comparison group
#    (rather than a fixed never-treated group)
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
  this.n <- nrow(gt_data)/3 # this is probably going to fail
  
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

  use_data <- this.data
  use_data$Y <- use_data$dY_base

  # after "adjusting" data, leverage calls to did
  attgt <- did_attgt(this.data, xformla, ...)
  
  attgt
}
