#' @title Setup data for covariate pre-test
#'
#' @description Sets up the data for conducting covariate pre-test
#'
#' @inheritParams setup_pte
#'
#' @return \code{pte_params} object
#'
#' @export
setup_cov_pre_test <- function(yname,
                               gname,
                               tname,
                               idname,
                               data,
                               nife=1,
                               anticipation=0,
                               ...) {

  required_pre_periods <- nife + 1
  
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
  #time.periods <- time.periods[1:(length(time.periods) - nife)]
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
                       cband=TRUE, # never use anything from here after
                       alp=0.05,
                       boot_type="multiplier",
                       biters=100,
                       cl=1)

  params
}
