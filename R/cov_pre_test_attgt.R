#' @title cov_pre_test_attgt
#' @description Computes estimates of ATT(g,t), overall ATT, and dynamic effects
#'  under interactive fixed effects model for untreated potential outcomes
#'  using the approach in Callaway and Karami (2021)
#'
#' @inheritParams ife_attgt
#' @param anticipation Number of periods that treatment is anticipated.  Default
#'  is 0.  This is in ``periods''; e.g., code will work in time periods are
#'  equally spaced but 2 years apart.  In this case, to allow for treatment
#'  anticipation of 2 year (<=> 1 period), set \code{anticipation = 1}.
#'
#' @return
cov_pre_test_attgt <- function(gt_data, nife=1,
                      xformla=~1, zformla, anticipation=0,  
                      ret_ife_regs=FALSE, ...) {

  if (nrow(gt_data) == 0) {
    return(attgt_if(attgt=NA, inf_func=NA, extra_gt_returns=NA))
  }
  # base period is the first one in this subset of the data
  base.period <- min(gt_data$period)
  tp <- max(gt_data$period)
  this.n <- nrow(gt_data)/(nife+2)
  
  # take difference with respect to base period
  this.data <- gt_data %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(dY_base=(Y-Y[period==base.period])) %>%
    as.data.frame()
  # and drop base period
  if (nife > 0) this.data <- subset(this.data, period != base.period)
  # add a lag variable to the data (hack!), just for the name
  this.data$.lag <- rev(sort(unique(gt_data$period)))[2]- this.data$period + 1
  
  # split pre and post data, eventually merge them back
  post.data <- subset(this.data, name == "post")
  post.data <- post.data %>% dplyr::rename(dY_post=dY_base)
  pre.data <- subset(this.data, name == "pre")

  # convert pre-data into cross-sectional data
  pre.data <- pre.data %>%
    select(id, .lag, dY_base) %>%
    dplyr::group_by(id) %>%
    tidyr::pivot_wider(names_prefix="dY_base",
                       names_from=.lag,
                       names_glue="dLagY{.lag}_base",
                       values_from=dY_base) %>%
    as.data.frame()

  # merge data, this is one row per unit and can use to run regressions
  # to identify ife model
  this.data <- dplyr::inner_join(post.data, pre.data, by="id")

  # hack to get extra column names for dY variables
  dY_names <- if (nife >  0) this.data %>% select(starts_with("dLagY")) %>% colnames else character(0)
  dY_names <- rev(dY_names)
  
  # formula for y ~ x
  outcome_formla <- BMisc::toformula(yname="dY_post", xnames=c(BMisc::rhs.vars(xformla), dY_names))
  zformla <- BMisc::toformula("", xnames=c(BMisc::rhs.vars(xformla), paste0("as.factor(G)")))
  
  # estimate ife model
  this.comparison <- subset(this.data, D==0)# subset(this.data, G != g)
  comparison_ids <- this.comparison$id
  comparison_p <- length(comparison_ids)/this.n
  ife_reg <- ivreg::ivreg(outcome_formla, instruments=zformla, data=this.comparison)
  
  
  extra_returns <- list(ife_reg=ife_reg)
  
  attgt_if(attgt=NULL, inf_func=NA, extra_gt_returns=extra_returns)
}
