#-----------------------------------------------------------------------------
# Some things to notice:
#
# * We rely on equally spaced periods
# * We use all not-yet-treated observations as the comparison group
#    (rather than a fixed never-treated group)
# * We difference out unobserved heterogeneity using (t-nife-1) as the
#    base period but other choices could work here
#-----------------------------------------------------------------------------


#' @title ife_attgt
#' @description Computes estimates of ATT(g,t), overall ATT, and dynamic effects
#'  under interactive fixed effects model for untreated potential outcomes
#'  using the approach in Callaway and Karami (2021)
#'
#' @param gt_data data frame that is local to the specific groups
#'  and times for which we'll be computing a treatment efect
#'  estimate.
#' @param nife The number of interactive fixed effects.  Default is 1.
#' @param xformla Formula for which covariates to include in the model.  Default is ~1.
#' @param zformla Formula for moment conditions to identify interactive fixed effects
#'  parameters.  This must include at least \code{nife} additional covariates relative to
#'  xformla.
#' @param ret_ife_regs Whether or not to return the first stage ife regressions; default is FALSE.
#' @param anticipation Number of periods that treatment is anticipated.  Default
#'  is 0.  This is in ``periods''; e.g., code will work in time periods are
#'  equally spaced but 2 years apart.  In this case, to allow for treatment
#'  anticipation of 2 year (<=> 1 period), set \code{anticipation = 1}.
#' @param ... extra arguments
#'
#' @return ptetools::attgt_if object
ife_attgt <- function(gt_data, nife=1,
                      xformla=~1, zformla, anticipation=0,  
                      ret_ife_regs=FALSE, ...) {

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
  
  # estimate ife model
  this.comparison <- subset(this.data, D==0)# subset(this.data, G != g)
  comparison_ids <- this.comparison$id
  comparison_p <- length(comparison_ids)/this.n
  ife_reg <- ivreg::ivreg(outcome_formla, instruments=zformla, data=this.comparison)
  #ife_reg <- estimatr::iv_robust
  # get the influence function from the first step
  first_step_if <- as.matrix(sandwich::estfun(ife_reg))
  first_step_if <- first_step_if %*% sandwich::bread(ife_reg)
  first_step_bet <- coef(ife_reg)

  #V <- bread(ife_reg) %*% (t(first_step_if) %*% first_step_if / ife_reg$n) %*% bread(ife_reg)
  # get attgt
  attgt_i <- subset(this.data, D==1)$dY_post - predict(ife_reg, newdata=subset(this.data, D==1))
  attgt <- mean(attgt_i)

  # get influence function for this part too
  this.treated <- subset(this.data, D==1)
  treated_ids <- this.treated$id
  treated_p <- length(treated_ids)/this.n
  mdY_post <- mean(this.treated$dY_post)
  this.treated_if1 <- this.treated$dY_post - mdY_post
  this.treated_Xmat <- model.matrix(outcome_formla, data=this.treated)
  mX <- apply(this.treated_Xmat, 2, mean)
  Xresid <- this.treated_Xmat - matrix(rep(mX, nrow(this.treated)), ncol=ncol(this.treated_Xmat), byrow=TRUE)
  this.treated_if2 <- as.matrix(Xresid) %*% as.matrix(first_step_bet)
  # account for not using the full sample
  second_step_if <- this.treated_if1 - this.treated_if2
  second_step_if <- second_step_if / treated_p

  # account for using estimating first step
  first_step_if <- -first_step_if %*% as.matrix(mX)
  # account for not using the full sample
  first_step_if <- first_step_if / comparison_p

  # put into influence function at the right spot
  this.if <- rep(0,this.n)
  idlist <- post.data$id
  this.if[idlist %in% comparison_ids] <- first_step_if
  this.if[idlist %in% treated_ids] <- second_step_if


  # model selection
  # cross-validation criteria
  Y <- ife_reg$y
  bet <- as.matrix(ife_reg$coefficients)
  X <- model.matrix(ife_reg$terms$regressors, data=ife_reg$model)
  Z <- model.matrix(ife_reg$terms$instruments, data=ife_reg$model)
  P <- Z %*% solve( t(Z) %*% Z) %*% t(Z)
  PX <- P %*% X
  PZ <- P %*% Z
  PP <- X %*% solve( t(PX) %*% PX) %*% t(PX)
  h <- diag(PP)
  ehat <- Y - X %*% bet
  eloo <- ehat / (1-h)
  cv_untreated <- mean(eloo^2)

  ## # K fold cross validation for treated group
  ## nu <- nrow(this.comparison)
  ## fold <- sample(1:5, size=nu, replace=TRUE)
  ## eloo <- rep(NA, nu)
  ## for (i in 1:5) {
  ##   kfold.comparison <- fold!=i
  ##   kfold <- fold==i
  ##   kfold.iv.coef <- coef(ivreg::ivreg(Y[kfold.comparison] ~ X[kfold.comparison,], ~Z[kfold.comparison,]))
  ##   kfold.iv.coef <- as.matrix(na.omit(kfold.iv.coef))
  ##   eloo[kfold] <- Y[kfold] - X[kfold,,drop=FALSE]%*%kfold.iv.coef
  ## }
  ## cv_untreated <- mean(eloo^2)
  
  # cross-validation for treated (this is useful in pre-treatment periods)
  cv_treated <- mean(attgt_i^2)

  #browser()
  
  # bayesian information criteria
  #Z <- model.matrix(zformla, data=this.comparison)
  u <- ife_reg$residuals
  n <- nrow(Z)
  k <- ncol(X)
  l <- ncol(Z)
  W <- (1/n) * t(Z) %*% Z
  gbar <- as.matrix(apply( (Z*u), 2, mean))
  J <-  n * t(gbar) %*% solve(W) %*% (gbar)

  bic <- J - log(n)*(l-k+1) # k already includes an extra term for each IFE #J - log(n)*(l - nife - k)#log(n)*(-nife)#0.75*nrow(W)#0.75( (q-p)(T-p) - k)

  if (!ret_ife_regs) {
    ife_reg <- NULL
  }
  
  extra_returns <- list(ife_reg=ife_reg, eloo=eloo, cv_untreated=cv_untreated, cv_treated=cv_treated, bic=bic)
  
  attgt_if(attgt, inf_func=this.if, extra_gt_returns=extra_returns)
}
