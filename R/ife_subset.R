#' @title ife_subset
#'
#' @description A function for obtaining the correct subset in a staggered treatment
#' adoption setup with panel data.  This is the subset of not yet treated
#' observations in a particular time period and in the previous
#' \code{nife+1} periods where \code{nife} is the number of interactive fixed effects.  
#'
#' @param data the full dataset
#' @param g the current group
#' @param tp the current time period
#' @param nife the number of interactive fixed effects
#' @param ... extra arguments to get the subset correct
#'
#' @return list that contains correct subset of data, \code{n1}
#'  number of observations
#'  in this subset, and \code{disidx} a vector of the correct ids for this
#'  subset.
#'
#' @export
ife_subset <- function(data, g, tp, nife, anticipation=0, ...) {

  # pre-treatment period used to difference out unobserved heterogeneity (w/o ife)
  main.base.period <- g - nife - 1 - anticipation
  

  #----------------------------------------------------
  # if it's a pre-treatment time period (used for the
  # pre-test, we need to adjust the base period)

  # group not treated yet
  if (tp < g) {
    # adjust base period earlier (relative to pre-treatment period)
    # don't no anticipation or anything here
    base.period <- tp - nife - 1
  } else {
    # this is a post-treatment period
    base.period <- main.base.period
  }
  #----------------------------------------------------

  # get group g and not-yet-treated group
  this.data <- subset(data, G==g | G>tp | G==0)

  # get current period and base period data
  this.data <- subset(this.data, period %in% c(tp, seq(base.period, length.out=(nife+1))))

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
