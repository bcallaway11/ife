#-----------------------------------------------------------------------------
#
# This is a copy of ife_subset except that we don't need to retain all
# intermediate periods.
#
#-----------------------------------------------------------------------------


#' @title lt_subset
#'
#' @description A function for obtaining the correct subset in a staggered treatment
#' adoption setup with panel data using a linear trends model for untreated potential
#' outcomes.
#'
#' @inheritParams ife_subset
#'
#' @return list that contains correct subset of data, \code{n1}
#'  number of observations
#'  in this subset, and \code{disidx} a vector of the correct ids for this
#'  subset.
#'
#' @export
lt_subset <- function(data, g, tp, anticipation=0, ...) {


  # pre-treatment period used to difference out unobserved heterogeneity (w/o ife)
  main.base.period <- g - 2 - anticipation
  

  #----------------------------------------------------
  # if it's a pre-treatment time period (used for the
  # pre-test, we need to adjust the base period)

  # group not treated yet
  if (tp < g) {
    # adjust base period earlier (relative to pre-treatment period)
    # don't no anticipation or anything here
    base.period <- tp - 2
  } else {
    # this is a post-treatment period
    base.period <- main.base.period
  }
  #----------------------------------------------------

  # get group g and not-yet-treated group
  this.data <- subset(data, G==g | G>tp | G==0)

  # get current period and base periods
  this.data <- subset(this.data, period %in% c(base.period, base.period+1, tp))

  # variable to keep track of pre/post periods
  this.data$name <- ifelse(this.data$period==tp, "post", "pre")

  # variable to indicate local treatment status
  this.data$D <- 1*(this.data$G==g)
  
  # make this.data into gt_data_frame object
  this.data <- gt_data_frame(this.data)

  # number of observations used for this (g,t)
  n1 <- length(unique(this.data$id))
  disidx <- unique(data$id) %in% unique(this.data$id)

  list(gt_data=this.data, n1=n1, disidx=disidx)
}
