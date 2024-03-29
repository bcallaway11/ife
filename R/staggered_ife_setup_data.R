#-----------------------------------------------------------------------------
#
# This is copy and paste of pte::setup_data, but adjusted to drop some
# periods at the end (where we do not have enough groups left to identify
# the model)
#
#-----------------------------------------------------------------------------

#' @title staggered_ife_setup_data
#'
#' @description Estimate treatment effects in an interactive fixed effects
#'  model for untreated potential outcomes by exploiting staggered treatment
#'  adoption.
#'
#' @inheritParams ife_setup_pte
#' @inheritParams pte::setup_pte
#' @param required_pre_periods the number of pre-treatment time
#'  periods that are needed
#' @param cband whether to compute a uniform (instead of pointwise)
#'  confidence band
#' @param boot_type type of bootstrap
#'
#' @return \code{pte_params} object
#' @export
staggered_ife_setup_pte <- function(yname,
                                    gname,
                                    tname,
                                    idname,
                                    data,
                                    required_pre_periods=1,
                                    nife=1,
                                    anticipation=0,
                                    cband=TRUE,
                                    alp=0.05,
                                    boot_type=boot_type,
                                    biters=100,
                                    cl=1,
                                    ...) {
  
  data <- as.data.frame(data)
  
  # setup data
  G <- data[,gname]
  id <- data[,idname]
  period <- data[,tname]

  data$G <- G
  data$id <- id
  n <- length(unique(data$id))
  # data$original_period <- period
  # data$original_group <- G
  data$Y <- data[,yname]
  
  time.periods <- unique(period)
  groups <- unique(data$G)

  # drop never treated group
  groups <- sort(groups)[-1]

  # do some recoding to make sure time periods are 1 unit apart
  # and then put these back together at the end
  # (this probably won't work for unequally spaced periods)
  original_time.periods <- time.periods
  original_groups <- groups

  # get new time periods / groups
  time.periods <- sapply(original_time.periods,
                         orig2t,
                         original_time.periods=original_time.periods)
  groups <- sapply(original_groups,
                   orig2t,
                   original_time.periods=original_time.periods)

  # put the period in new time scale
  data$period <- sapply(period,
                        orig2t,
                        original_time.periods=original_time.periods)
  data$G <- sapply(G,
                   orig2t,
                   original_time.periods=original_time.periods)

  # sort the time periods and drop the first
  # \code{required_pre_periods} time periods
  # these are the ones we loop over below
  time.periods <- sort(time.periods)[-seq(1,required_pre_periods)]

  #-----------------------------------------------------------------------------
  # this is what's different relative to `setup_pte`
  #-----------------------------------------------------------------------------
  last_attgt_period_idx <- tail(which(sapply(time.periods, function(tp) sum(tp<groups)) >= nife),1)
  time.periods <- time.periods[1:last_attgt_period_idx]
  #-----------------------------------------------------------------------------
  
  groups <- groups[groups %in% time.periods]
  # account for anticipation
  groups <- groups[ groups >= (min(time.periods)+anticipation) ]

  # drop early treated groups
  data <- data[ ! (data$G %in% seq(1,required_pre_periods+anticipation)), ]

  params <- pte_params(yname=yname,
                       gname=gname,
                       tname=tname,
                       idname=idname,
                       data=data,
                       glist=groups,
                       tlist=time.periods,
                       cband=cband,
                       alp=alp,
                       boot_type=boot_type,
                       biters=biters,
                       cl=cl)

  params
}

#' @title staggered_ife_setup_data2
#'
#' @description This function sets up the periods/groups for using 
#'  `staggered_ife_attgt2` (the main approach proposed in Callaway and 
#'  Tsyawo (2023)).
#'
#' @inheritParams ife_setup_pte
#' @inheritParams pte::setup_pte
#' @param required_pre_periods the number of pre-treatment time
#'  periods that are needed
#' @param cband whether to compute a uniform (instead of pointwise)
#'  confidence band
#' @param boot_type type of bootstrap
#'
#' @return \code{pte_params} object
#' @export
staggered_ife_setup_data2 <- function(yname,
                                    gname,
                                    tname,
                                    idname,
                                    data,
                                    required_pre_periods=1,
                                    nife=1,
                                    anticipation=0,
                                    cband=TRUE,
                                    alp=0.05,
                                    boot_type=boot_type,
                                    biters=100,
                                    cl=1,
                                    ...) {
  
  data <- as.data.frame(data)
  
  # setup data
  G <- data[,gname]
  id <- data[,idname]
  period <- data[,tname]
  
  data$G <- G
  data$id <- id
  n <- length(unique(data$id))
  # data$original_period <- period
  # data$original_group <- G
  data$Y <- data[,yname]
  
  time.periods <- unique(period)
  groups <- unique(data$G)
  
  # drop never treated group
  groups <- sort(groups)[-1]
  
  # do some recoding to make sure time periods are 1 unit apart
  # and then put these back together at the end
  # (this probably won't work for unequally spaced periods)
  original_time.periods <- time.periods
  original_groups <- groups
  
  # get new time periods / groups
  time.periods <- sapply(original_time.periods,
                         orig2t,
                         original_time.periods=original_time.periods)
  groups <- sapply(original_groups,
                   orig2t,
                   original_time.periods=original_time.periods)
  
  # put the period in new time scale
  data$period <- sapply(period,
                        orig2t,
                        original_time.periods=original_time.periods)
  data$G <- sapply(G,
                   orig2t,
                   original_time.periods=original_time.periods)
  
  # sort the time periods and drop the first
  # \code{required_pre_periods} time periods
  # these are the ones we loop over below
  time.periods <- sort(time.periods)[-seq(1,required_pre_periods)]
  
  # include first period if nife==-1 (=> level comparisons)
  if (nife==-1) {
    time.periods <- c(1,time.periods)
  }
  
  #-----------------------------------------------------------------------------
  # this is what's different relative to `setup_pte`
  #-----------------------------------------------------------------------------
  last_attgt_period_idx <- tail(which(sapply(time.periods, function(tp) sum(tp<groups)) >= nife),1)
  time.periods <- time.periods[1:last_attgt_period_idx]
  #-----------------------------------------------------------------------------
  
  groups <- groups[groups %in% time.periods]
  # account for anticipation
  groups <- groups[ groups >= (min(time.periods)+anticipation) ]
  
  
  #----------------------------------------------------------------
  #  this is the difference relative to `staggered_ife_setup_data`,
  #  we are not going to drop any data
  #----------------------------------------------------------------
  #  data <- data[ ! (data$G %in% seq(1,required_pre_periods+anticipation)), ]
  
  params <- pte_params(yname=yname,
                       gname=gname,
                       tname=tname,
                       idname=idname,
                       data=data,
                       glist=groups,
                       tlist=time.periods,
                       cband=cband,
                       alp=alp,
                       boot_type=boot_type,
                       biters=biters,
                       cl=cl)
  
  params
}