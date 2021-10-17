
#' @title ife_setup_pte
#'
#' @description Setup panel treatment effects parameters
#'
#' @inheritParams pte::pte_params
#' @inheritParams ife
#' @param ... extra arguments
#'
#' @export
ife_setup_pte <- function(yname,
                          gname,
                          tname,
                          idname,
                          data,
                          nife,
                          anticipation=0,
                          alp=0.05,
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
  data$period <- period
  data$Y <- data[,yname]
  
  time.periods <- unique(data$period)
  groups <- unique(data$G)

  # drop never treated group
  groups <- sort(groups)[-1]

  # do some recoding to make sure time periods are 1 unit apart
  # and then put these back together at the end
  # (this probably won't work for unequally spaced periods)
  original_time.periods <- time.periods
  original_groups <- groups

  # new time periods
  unique_t <- seq(1,length(unique(time.periods)))
  # function to switch from "new" t values to  original t values
  t2orig <- function(t) {
    unique(c(original_time.periods,0))[which(c(unique_t,0)==t)]
  }
  # function to switch between "original"
  #  t values and new t values
  orig2t <- function(orig) {
    c(unique_t,0)[which(unique(c(original_time.periods,0))==orig)]
  }
  # get new time periods / groups
  time.periods <- sapply(original_time.periods, orig2t)
  groups <- sapply(original_groups, orig2t)

  
  # put the period in new time scale
  data$period <- sapply(data$period, orig2t)
  data$G <- sapply(data$G, orig2t)

  # sort the time periods and drop the first (nife+1) time periods
  # these are the ones we loop over below
  time.periods <- sort(time.periods)[-seq(1,nife+1)]
  groups <- groups[groups %in% time.periods]
  # account for anticipation
  groups <- groups[ groups >= (min(time.periods)+anticipation) ] 

  params <- pte_params(yname=yname,
                       gname=gname,
                       tname=tname,
                       idname=idname,
                       data=data,
                       glist=groups,
                       tlist=time.periods,
                       alp=alp,
                       biters=biters,
                       cl=cl)

  params
}
