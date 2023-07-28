#------------------------------------------------------------------------
#  These are functions used for testing in the `ife` package.
#  Many of the default choices come from the simulations in 
#  Callaway and Karami (2023) and Callaway and Tsyawo (2023).
#------------------------------------------------------------------------

#' @title reset.sim A function that sets "reasonable" values for parameters 
#'  with which to simulate data using the `gen_ife_data` function
#' @keywords internal
#' @export
reset.sim2 <- function() {
  nobs <- 1000
  F1 <- c(1,-log(1),1^2)
  F2 <- c(2,2*log(2),2^2)
  F3 <- c(3,-3*log(3),3^2)
  F4 <- c(4,4*log(4),4^2)
  F5 <- c(5,-5*log(5),-3^2)
  F6 <- c(6,6*log(6),-2^2)
  F7 <- c(7,-7*log(7),-1^2) 
  nife <- 1
  nife_actual <- 1
  rho <- c(1,1,1)
  
  list(nobs=nobs, 
       F1=F1, 
       F2=F2, 
       F3=F3, 
       F4=F4, 
       F5=F5, 
       F6=F6, 
       F7=F7, 
       nife=nife, 
       nife_actual=nife_actual, 
       rho=rho)
}


#' @keywords internal
#' 
#' @export
gen_ife_data2 <- function(sim_params_list) {
  nobs <- sim_params_list$nobs
  F1 <- sim_params_list$F1
  F2 <- sim_params_list$F2
  F3 <- sim_params_list$F3
  F4 <- sim_params_list$F4
  F5 <- sim_params_list$F5
  F6 <- sim_params_list$F6
  F7 <- sim_params_list$F7
  nife <- sim_params_list$nife
  nife_actual <- sim_params_list$nife_actual
  rho <- sim_params_list$rho
  
  n <- nobs
  p <- rep(1/5, 4) # repeat 1/5 4 times
  
  nZ <- max(nife, 0)
  
  # there are three periods before any unit treated, and then groups from 4 to 7 + never-treated
  G <- sample(seq(5,8), size=n, replace=TRUE, prob=p)
  G[G==8] <- 0 # "never treated"
  
  # draw eta and lambda to be correlated with group
  eta <- G + rnorm(n, sd=sqrt(.1))
  if (nife_actual == -1) eta <- eta - G # unit fe not correlated with group in this case
  Z1 <- rnorm(n)
  Z2 <- rnorm(n)
  Z3 <- rnorm(n)
  lam1 <- 1 + 2*G + rho[1]*Z1 + rnorm(n, sd=sqrt(.1))
  lam2 <- 1 - 5*G + rho[2]*Z2 + rnorm(n, sd=sqrt(.1))
  lam3 <- 5 - 10*G + rho[3]*Z3 + rnorm(n, sd=sqrt(.1))
  #
  
  F2 <- F2
  F3 <- F3
  U1 <- rnorm(n, sd=sqrt(.1))#rt(n,5) #rnormmix(n, rep(.5,2), mu=c(-2,2), sigma=1)
  U2 <- rnorm(n, sd=sqrt(.1))#rt(n,5)#rnormmix(n, rep(.5,2), mu=c(-2,2), sigma=1)
  U3 <- rnorm(n, sd=sqrt(.1))#rt(n,5)#rnorm(n)
  U4 <- rnorm(n, sd=sqrt(.1))#rt(n,5)
  U5 <- rnorm(n, sd=sqrt(.1))#rt(n,5)
  U6 <- rnorm(n, sd=sqrt(.1))
  U7 <- rnorm(n, sd=sqrt(.1))
  thet1 <- 0
  thet2 <- .1
  thet3 <- .2
  thet4 <- .3
  thet5 <- .4
  thet6 <- .5
  thet7 <- .6
  X1 <- 1*rnorm(n)
  X2 <- 2*rnorm(n)
  X3 <- 3*rnorm(n)
  X <- X1+X2  # X1+X2+X3

  Y01 <- thet1 + eta + X + U1
  Y02 <- thet2 + eta + X + U2
  Y03 <- thet3 + eta + X + U3
  Y04 <- thet4 + eta + X + U4
  Y05 <- thet5 + eta + X + U5
  Y06 <- thet6 + eta + X + U6
  Y07 <- thet7 + eta + X + U7
  if (nife_actual > 0) {
    Y01 <- Y01 + F1[1]*lam1
    Y02 <- Y02 + F2[1]*lam1
    Y03 <- Y03 + F3[1]*lam1
    Y04 <- Y04 + F4[1]*lam1
    Y05 <- Y05 + F5[1]*lam1
    Y06 <- Y06 + F6[1]*lam1
    Y07 <- Y07 + F7[1]*lam1
  }
  if (nife_actual > 1) {
    Y01 <- Y01 + F1[2]*lam2
    Y02 <- Y02 + F2[2]*lam2
    Y03 <- Y03 + F3[2]*lam2
    Y04 <- Y04 + F4[2]*lam2
    Y05 <- Y05 + F5[2]*lam2
    Y06 <- Y06 + F6[2]*lam2
    Y07 <- Y07 + F7[2]*lam2
  }
  if (nife_actual > 2) {
    Y01 <- Y01 + F1[3]*lam3
    Y02 <- Y02 + F2[3]*lam3
    Y03 <- Y03 + F3[3]*lam3
    Y04 <- Y04 + F4[3]*lam3
    Y05 <- Y05 + F5[3]*lam3
    Y06 <- Y06 + F6[3]*lam3
    Y07 <- Y07 + F7[3]*lam3
  }  
  
  Y1 <- Y01
  Y2 <- Y02
  Y3 <- Y03
  Y4 <- (G<=4)*Y04 + (G>4)*Y04 # no effect of treatment
  Y5 <- (G<=5)*Y05 + (G>5)*Y05 # no effect of treatment
  Y6 <- (G<=6)*Y06 + (G>6)*Y06 # no effect of treatment
  Y7 <- (G<=7)*Y07 + (G>7)*Y07 # no effect of treatment
  id <- seq(1:n)
  dta <- cbind.data.frame(G, 
                          Y1,Y2,Y3,Y4,Y5,Y6,Y7, 
                          X1,X2,
                          Z1,Z2,Z3, 
                          id)
  dta.expanded <- cbind.data.frame(G,
                                   Y1,Y2,Y3,Y4,Y5,Y6,Y7,
                                   eta,
                                   lam1,lam2,lam3,
                                   Z1,Z2,Z3,
                                   id)
  
  # convert to long/panel data format
  dta.long <- pivot_longer(dta, 
                           cols=starts_with("Y"), 
                           names_to="period", 
                           names_prefix="Y", 
                           names_transform=as.integer, 
                           values_to="Y") %>% as.data.frame()
  
  dta.long
}