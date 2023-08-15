#-----------------------------------------------------------------------------
# Some things to notice:
#
# * We rely on equally spaced periods
# * We use all not-yet-treated observations as the comparison group
#    (rather than a fixed never-treated group)
# * We difference out unobserved heterogeneity using (t-nife-1) as the
#    base period but other choices could work here
#-----------------------------------------------------------------------------


#' @title staggered_ife_attgt
#' @description Computes estimates of group-time average treatment
#'  effects in an interactive treatment effects model for untreated
#'  potential outcomes by exploiting staggered treatment adoption as 
#'  in Callaway and Tsyawo (2023).  This function is based on the local-
#'  differencing approach (similar to the approach proposed in 
#'  Callaway and Karami).  See `staggered_ife_attgt2` for the main approach
#'  discussed in the paper where all pre-treatment periods are used 
#'  to estimate the interactive fixed effects model.
#'
#' @inheritParams ife_attgt
#' @return \code{pte::attgt_if} object
#' @export
staggered_ife_attgt <- function(gt_data,
                                nife=1,
                                xformla=~1,
                                anticipation=0,  
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
  this.data <- subset(this.data, period != base.period)

  
  # split pre and post data, eventually merge them back
  post.data <- subset(this.data, name == "post")
  post.data <- post.data %>% dplyr::rename(dY_post=dY_base)

  if (nife > 0) {
    pre.data <- subset(this.data, name == "pre")
    
    # convert pre-data into cross-sectional data
    pre.data <- pre.data %>%
      select(id, period, dY_base) %>%
      dplyr::group_by(id) %>%
      tidyr::pivot_wider(names_prefix="dY_base", names_from=period, values_from=dY_base) %>%
      as.data.frame()
    
    # merge data, this is one row per unit and can use to run regressions
    # to identify ife model
    this.data <- dplyr::inner_join(post.data, pre.data, by="id")
  } else {
    this.data <- post.data
  }

  # hack to get extra column names for dY variables
  dY_names <- this.data %>% select(starts_with("dY_base")) %>% colnames
  
  # formula for y ~ x
  outcome_formla <- BMisc::toformula(yname="dY_post", xnames=c(BMisc::rhs.vars(xformla), dY_names))

  #-----------------------------------------------------------------------------
  # this is only change relative to ife_attgt
  if (nife > 0) {
    zformla <- BMisc::toformula(yname="", xnames=c(BMisc::rhs.vars(xformla), "as.factor(G)"))
  } else {
    zformla <- ~1
  }
  #-----------------------------------------------------------------------------

  
  # estimate ife model
  this.comparison <- subset(this.data, D==0)# subset(this.data, G != g)
  comparison_ids <- this.comparison$id
  comparison_p <- length(comparison_ids)/this.n

  ife_reg <- AER::ivreg(outcome_formla, instruments=zformla, data=this.comparison)
  # get the influence function from the first step
  first_step_if <- as.matrix(sandwich::estfun(ife_reg))
  first_step_if <- first_step_if %*% sandwich::bread(ife_reg)
  first_step_bet <- coef(ife_reg)

  #V <- bread(ife_reg) %*% (t(first_step_if) %*% first_step_if / ife_reg$n) %*% bread(ife_reg)
  # get attgt

  #-----------------------------------------------------------------------------
  # this is also different, due to having to hack AER
  #-----------------------------------------------------------------------------
  this.treated <- subset(this.data, D==1)
  this.treated$G <- unique(this.comparison$G)[1] # doesn't do anything, but stops AER from crashing
  attgt <- mean(subset(this.data, D==1)$dY_post) - mean(predict(ife_reg, newdata=this.treated))

  # get influence function for this part too
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

  if (!ret_ife_regs) {
    ife_reg <- NULL
  }
  
  attgt_if(attgt, inf_func=this.if, extra_gt_returns=ife_reg)
}



#' @title staggered_ife_attgt2
#' @description Computes estimates of group-time average treatment
#'  effects in an interactive treatment effects model for untreated
#'  potential outcomes by exploiting staggered treatment adoption as 
#'  in Callaway and Tsyawo (2023).  This function uses the main approach
#'  discussed in the paper where all pre-treatment periods are used 
#'  to estimate the interactive fixed effects model.
#'
#' @inheritParams ife_attgt
#' @inheritParams staggered_ife2
#' @return \code{pte::attgt_if} object
#' @export
staggered_ife_attgt2 <- function(gt_data,
                                nife=1,
                                weighting_matrix="gmm",
                                xformla=~1,
                                anticipation=0,  
                                ret_ife_regs=FALSE, ...) {
  
  tp <- max(gt_data$period)
  this.n <- length(unique(gt_data$id))
  this.g <- unique(subset(gt_data, D==1)$G)
  
  if (anticipation != 0) stop("anticipation is not yet supported")
  if ( !(xformla==~1) ) stop("including covariates is not yet supported")
  
  # figure out the base period...(g-1)
  base.period <- unique(subset(gt_data,D==1)$G)-1
  # if this is a pre-treatment period, pick period
  # right before it
  base.period <- min(base.period, tp-1)
  
  # list of available comparison groups
  gcomplist <- sort(unique(subset(gt_data, G!=this.g)$G))
  
  this.data <- gt_data
  this.data$dY <- BMisc::get_first_difference(gt_data, "id", "Y", "period")
  
  # handle case with nife==-1 (level-comparisons)
  if (nife == -1) {
    this.data$dY_base <- this.data$Y # this is hack to get level comparison
  } else {
    # main case!
    # take difference with respect to base period
    this.data$Ygmin1 <- get_Yit(this.data, 
                                tp=base.period,
                                idname="id", 
                                yname="Y", 
                                tname="period")
    this.data$dY_base <- this.data$Y - this.data$Ygmin1
    
    # drop first period due to first differencing
    this.data <- subset(this.data, period != min(this.data$period))
  }
    
  # convert to cross sectional data
  # split pre and post data, eventually merge them back
  post.data <- subset(this.data, name == "post")
  
  
  # main case, nife > 0
  if (nife > 0) {
    pre.data <- subset(this.data, name == "pre")
    pre.data <- pre.data %>%
      select(id, G, period, dY) %>%
      dplyr::group_by(id) %>%
      tidyr::pivot_wider(names_prefix="dY", names_from=period, values_from=dY) %>%
      as.data.frame()
    
    # get principal components jointly for all groups
    pre.data_pc <- prcomp(select(pre.data, starts_with("dY")), rank.=nife)$x
    pre.data <- cbind.data.frame(pre.data, pre.data_pc)
    
    pre.data_untreated <- subset(pre.data, G!=this.g)
    pre_untreated_pca_nife <- select(pre.data_untreated, starts_with("PC")) %>% as.matrix()
    # add intercept
    pre_untreated_pca_nife <- cbind(1, pre_untreated_pca_nife)
    
    #pre_untreated_pca_nife <- prcomp(select(pre.data_untreated, starts_with("dY")),
    #                                 rank.=nife)$x
    #pre.data_untreated <- cbind.data.frame(pre.data_untreated, pre_untreated_pca_nife)
    
    Gamma_gt <- pre.data_untreated %>%
      select(G, starts_with("PC")) %>%
      group_by(G) %>%
      summarize_all("mean") %>% 
      select(starts_with("PC")) %>% as.matrix()
    
    #Gamma_gt <- pre.data_untreated %>%
    #  select(G, starts_with("dY")) %>%
    #  group_by(G) %>% 
    #  summarize(across(starts_with("dY"), mean)) %>%
    #  select(starts_with("dY")) %>% as.matrix()
    
    Gamma_gt <- cbind(1, Gamma_gt)
  }
  
  post.data_untreated <- subset(post.data, D==0)
  n_untreated <- nrow(post.data_untreated)
  LdY_base <- post.data_untreated %>% 
    group_by(G) %>%
    dplyr::summarize(dY_base=mean(dY_base))
  LdY_base <- LdY_base$dY_base
  
  # handle case w nife==0
  if (nife %in% c(0,-1)) {
    # this creates a vector of 1's 
    # I think this is ok, but it is different from above,
    # earlier it was convenient not to divide by pg; here it is convenient
    # to divide by pg which turns every element to be a 1.
    Gamma_gt <- matrix(rep(1,length(gcomplist)), ncol=1)
    pre_untreated_pca_nife <- matrix(rep(1,n_untreated), ncol=1)# name is awkward, but the only regressor here is intercept
  }

  #------------------------------------------------------------------------
  #  first-step GMM estimation
  #------------------------------------------------------------------------
  
  # handle bug that happens in nife==0 case for last-treated group
  if (length(unique(post.data_untreated$G))==1) {
    Z_untreated <- matrix(rep(1,n_untreated), ncol=1)
  } else {
    # this is the main case
    Z_untreated <- model.matrix(~-1 + as.factor(G), data=post.data_untreated)
  }
  
  # settle on weighting matrix
  if (weighting_matrix == "identity") {
    W <- diag(nrow=ncol(Z_untreated))
  } else { # 2sls weighting matrix
    W <- solve(t(Z_untreated)%*%Z_untreated/n_untreated) # 2sls weighting matrix
  }
  
  first_step_params <- solve( t(Gamma_gt) %*% W %*% Gamma_gt ) %*% t(Gamma_gt) %*% W %*% LdY_base
  first_step_yhat <- pre_untreated_pca_nife %*% first_step_params
  first_step_ehat <- post.data_untreated$dY_base - first_step_yhat
  
  Ze_untreated <- Z_untreated*as.numeric(first_step_ehat)
  
  # two-step gmm, if requested
  if (weighting_matrix == "gmm") {
    W <- t(Ze_untreated)%*%Ze_untreated/n_untreated
    
    # same code, new weighting matrix
    first_step_params <- solve( t(Gamma_gt) %*% W %*% Gamma_gt ) %*% t(Gamma_gt) %*% W %*% LdY_base
    first_step_yhat <- pre_untreated_pca_nife %*% first_step_params
    first_step_ehat <- post.data_untreated$dY_base - first_step_yhat
  }
  
  # calculate influence function for first-step estimation
  first_step_if <- t(solve( t(Gamma_gt) %*% W %*% Gamma_gt ) %*% t(Gamma_gt) %*% W %*% t(Ze_untreated))
  
  # if requested, we'll return the first-step estimates, code to do it:
  first_step_se <- sqrt(diag(t(first_step_if)%*%first_step_if)) / sqrt(n_untreated)
  
  ife_reg <- list(coefs=first_step_params, 
                  se=first_step_se)
  
  #------------------------------------------------------------------------
  #  estimate att(g,t) given the first-step estimates of parameters
  #  from the interactive fixed effects model
  #------------------------------------------------------------------------
  post.data_treated <- subset(post.data, D==1)
  n_treated <- nrow(post.data_treated)
  
  # main case nife > 0
  if (nife > 0) {
    pre.data_treated <- subset(pre.data, G==this.g)
    pca_nife <- select(pre.data_treated, starts_with("PC")) %>% as.matrix()
  #pca_nife <- prcomp(select(pre.data_treated, starts_with("dY")),
  #       rank.=nife)$x
    pca_nife <- cbind(1, pca_nife)
  }
  
  if (nife %in% c(0,-1)) pca_nife <- matrix(rep(1,n_treated))
  
  this.attgt <- mean(post.data_treated$dY_base) - mean(pca_nife %*% first_step_params)
  
  # compute influence function for second step
  pg <- mean(this.data$G==this.g)
  second_step_if1 <- post.data_treated$dY_base - mean(post.data_treated$dY_base)
  m_pca_nife <- apply(pca_nife, 2, mean)
  pca_resid <- pca_nife - matrix(rep(m_pca_nife, n_treated), ncol=length(m_pca_nife), byrow=TRUE)
  second_step_if2 <- pca_resid %*% first_step_params
  second_step_if <- second_step_if1 - second_step_if2
  second_step_if <- second_step_if / pg # accounts for only using treated group in this step
  
  # adjust first step influence function to account for where it enters
  # expression on ATT(g,t)
  first_step_if <- -(first_step_if %*% as.matrix(m_pca_nife)) / (1-pg) 
  
  #browser()
  
  # estimating pg component of variance
  # pg_if <- this.attgt*( (1*(post.data$G==this.g)) - pg) / pg
  
  # set up the overall influence function to return
  inf_func <- rep(0, this.n)
  idlist <- post.data$id
  treated_ids <- post.data_treated$id
  comparison_ids <- post.data_untreated$id
  inf_func[idlist %in% comparison_ids] <- as.numeric(first_step_if)
  inf_func[idlist %in% treated_ids] <- as.numeric(second_step_if)
  # inf_func <- inf_func - pg_if
    
  #browser()
  
  #------------------------------------------------------------------------
  #  check standard errors manually, can delete this later
  #------------------------------------------------------------------------
  V <- t(inf_func) %*% inf_func / this.n
  se <- sqrt(V)/sqrt(this.n)
  se
  V1 <- t(first_step_if)%*%first_step_if / this.n
  se1 <- sqrt(V1)/sqrt(this.n)
  V2 <- t(second_step_if)%*%second_step_if / this.n
  se2 <- sqrt(V2)/sqrt(this.n)
  se1
  se2
  V1 + V2
  
  # #------------------------------------------------------------------------
  # #  some model selection code from ife function...not currently used
  # #------------------------------------------------------------------------
  # # model selection
  # # cross-validation criteria
  # 
  # Y <- ife_reg$y
  # bet <- as.matrix(ife_reg$coefficients)
  # X <- model.matrix(ife_reg$terms$regressors, data=ife_reg$model)
  # Z <- model.matrix(ife_reg$terms$instruments, data=ife_reg$model)
  # P <- Z %*% solve( t(Z) %*% Z) %*% t(Z)
  # PX <- P %*% X
  # PZ <- P %*% Z
  # PP <- X %*% solve( t(PX) %*% PX) %*% t(PX)
  # h <- diag(PP)
  # ehat <- Y - X %*% bet
  # eloo <- ehat / (1-h)
  # cv_untreated <- mean(eloo^2)
  # 
  # ## # K fold cross validation for treated group
  # ## nu <- nrow(this.comparison)
  # ## fold <- sample(1:5, size=nu, replace=TRUE)
  # ## eloo <- rep(NA, nu)
  # ## for (i in 1:5) {
  # ##   kfold.comparison <- fold!=i
  # ##   kfold <- fold==i
  # ##   kfold.iv.coef <- coef(ivreg::ivreg(Y[kfold.comparison] ~ X[kfold.comparison,], ~Z[kfold.comparison,]))
  # ##   kfold.iv.coef <- as.matrix(na.omit(kfold.iv.coef))
  # ##   eloo[kfold] <- Y[kfold] - X[kfold,,drop=FALSE]%*%kfold.iv.coef
  # ## }
  # ## cv_untreated <- mean(eloo^2)
  # 
  # # cross-validation for treated (this is useful in pre-treatment periods)
  # cv_treated <- mean(attgt_i^2)
  # 
  # #browser()
  # 
  # # bayesian information criteria
  # #Z <- model.matrix(zformla, data=this.comparison)
  # u <- ife_reg$residuals
  # n <- nrow(Z)
  # k <- ncol(X)
  # l <- ncol(Z)
  # W <- (1/n) * t(Z) %*% Z
  # gbar <- as.matrix(apply( (Z*u), 2, mean))
  # J <-  n * t(gbar) %*% solve(W) %*% (gbar)
  # 
  # bic <- J - log(n)*(l-k+1) # k already includes an extra term for each IFE #J - log(n)*(l - nife - k)#log(n)*(-nife)#0.75*nrow(W)#0.75( (q-p)(T-p) - k)
  # 
  # if (!ret_ife_regs) {
  #   ife_reg <- NULL
  # }
  # 
  # extra_returns <- list(ife_reg=ife_reg, eloo=eloo, cv_untreated=cv_untreated, cv_treated=cv_treated, bic=bic)
  
  extra_returns <- list(ife_reg=ife_reg)
  
  attgt_if(this.attgt, inf_func=inf_func, extra_gt_returns=extra_returns)
}
