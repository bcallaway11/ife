#' @title compute.ife
#' @description Does the heavy-lifting of actually computing treatment effects
#'
#' @param formla y ~ d where y is the outcome and d is the treatment indicator
#' @param xformla ~x where x includes the covariates in the model
#' @param data is the dataset
#' @param treat.period the numerical value of the time period when individuals become treated
#'
#' @return IFEobj
#' @keywords internal
#' @export
compute.ife <- function(formla, xformla, data, idname, tname, treat.t) {
    yname <- BMisc::lhs.vars(formla)
    dname <- BMisc::rhs.vars(formla)
    data$y <- data[,BMisc::lhs.vars(formla)]
    data$d <- data[,BMisc::rhs.vars(formla)]
    ##figure out the dates and make balanced panel
    tlist <- unique(data[,tname])[order(unique(data[,tname]))] ## this is going to be from smallest to largest

    ##flist <- unique(data[,first.treat.name])[order(unique(data[,first.treat.name]))]
    ##flist <- flist[flist>0]

    ##################################
    ## do some error checking
    if (!is.numeric(tlist)) {
        warning("not guaranteed to order time periods correclty if they are not numeric")
    }

    if (!(length(dname)==1)) {
        stop("you should only include one treatment variable")
    }

    ## check that all treatments are either 0 or 1
    ####################################

    
    tlen <- length(tlist)

    ####################################
    ## collect all the things that don't change across each period
    Y1U <- data[data[,tname]==tlist[1] & data[,dname]==0, "y"]
    Y2U <- data[data[,tname]==tlist[2] & data[,dname]==0, "y"]
    Y3U <- data[data[,tname]==tlist[3] & data[,dname]==0, "y"]
    dta1U <- data[data[,tname]==tlist[1] & data[,dname]==0, ] ## this is ruling out time varying regressors
    XU <- model.matrix(xformla, dta1U)
    Y1T <- data[data[,tname]==tlist[1] & data[,dname]==1, "y"]
    Y2T <- data[data[,tname]==tlist[2] & data[,dname]==1, "y"]
    Y3T <- data[data[,tname]==tlist[3] & data[,dname]==1, "y"]
    dta1T <- data[data[,tname]==tlist[1] & data[,dname]==1, ]
    XT <- model.matrix(xformla, dta1T)

    
    ## get some parts of influence function
    XTbar <- apply(XT, 2, mean)
    ifXT <- t(apply(XT, 1, function(x) x  - XTbar))
    ifY12T <- Y2T - Y1T - mean(Y2T-Y1T)
    Y12Tbar <- mean(Y2T-Y1T)

    ATTout <- c() ## get results for ATT
    IFout <- c() ## get results for the influence function
    for (i in 4:tlen) {
        ## get 1st step iv results
        YtU <- data[data[,tname]==tlist[i] & data[,dname]==0, "y"]
        ZZU <- cbind(Y3U,XU)
        XXU <- cbind(Y2U - Y1U, XU)
        YYU <- as.matrix(YtU - Y1U)
        bet <- solve(t(ZZU)%*%XXU)%*%t(ZZU)%*%YYU
        delt <- bet[-1]
        Ft <- bet[1]

        ## calculate ATT
        YtT <- data[data[,tname]==tlist[i] & data[,dname]==1, "y"]
        XTbar <- as.matrix(apply(XT, 2, mean))
        ATT <- mean(YtT - Y1T) - t(XTbar)%*%delt - Ft*mean(Y2T - Y1T)

        ATTout[length(ATTout)+1] <- ATT

        ##browser()
        ## calculate influence function
        ## if1 <- (YtT - Y1T) - mean(YtT - Y1T)
        ## ifbet <-
        ##     solve(t(ZZU)%*%XXU)
        ##     as.matrix(ZZU*as.numeric(YYU - XXU%*%bet))
        ## if2 <- t(XTbar)
        
    }

    ###########################################################################
    ## for time period 3, just use time period 4 as instrument in 1st step
    ## handle it as a special case here and put it in the first place in ATTout
    Y4U <- data[data[,tname]==tlist[4] & data[,dname]==0, "y"]
    ZZU <- cbind(Y4U,XU)
    XXU <- cbind(Y2U - Y1U, XU)
    YYU <- as.matrix(Y3U - Y1U)
    bet <- solve(t(ZZU)%*%XXU)%*%t(ZZU)%*%YYU
    delt <- bet[-1]
    Ft <- bet[1]
    
    ## calculate ATT
    Y3T <- data[data[,tname]==tlist[3] & data[,dname]==1, "y"]
    XTbar <- as.matrix(apply(XT, 2, mean))
    ATT <- mean(Y3T - Y1T) - t(XTbar)%*%delt - Ft*mean(Y2T - Y1T)
    ###########################################################################
    

    ATTout <- c(ATT, ATTout)

    ##########################################################################
    ## next, we need to calculate standard errors
    ##if1  <- 

    #########################################################################
    ATTout
}


#' @title ife
#'
#' @description compute treatment effects in interactive fixed effects models
#' @inheritParams compute.ife
#'
#' @return IFEobj
#' @export
ife  <- function(formla, xformla, data, idname, tname, treat.t, boot.iters=100,
                 cores=1) {
    cat("Computing point estimates...")
    cat("Done\n")
    res <- compute.ife(formla, xformla, data, idname, tname, treat.t)
    cat("Bootstrapping standard errors...\n")
    boot.out <- pbapply::pbsapply(1:boot.iters, function(i) {
        bdta <- BMisc::blockBootSample(data, idname)
        compute.ife(formla, xformla, bdta, idname, tname, treat.t)
    }, cl=cores)
    boot.out <- t(boot.out)
    V <- var(boot.out) ## t(boot.out)%*%boot.out/boot.iters
    ##se <- apply(boot.out, 2, sd)
    n <- length(unique(data[,idname]))
    se <- diag(V)/sqrt(n)
    IFEobj(att=res, se=se, V=V)
}

#' @title IFEobj
#'
#' @description object that holds results of treatment effects in interactive
#'  fixed effects models
#'
#' @return IFEobj
#' @export
IFEobj <- function(att, se=NULL, V=NULL) {
    out <- list()
    out$att <- att
    out$se <- se
    out$V <- V
    class(out) <- "IFEobj"
    out
}


#' @title ggife
#'
#' @description plots treatment effects in interactive fixed effects models
#'
#' @param ifeobj an interactive fixed effects object
#'
#' @export
ggife <- function(ifeobj) {
    ## Idea: set the class of the ifeobj to be "did" so that we can leverage
    ## plotting capabilities in did package
    class(ifeobj) <- "did"
    ggdid(ifeobj)
}

#' @title compute.ife2
#'
#' @description does the heavy-lifting of estimating average treatment effects
#'  in interactive fixed effects models using GMM
#' @inheritParams compute.ife
#'
#' @return IFEobj
#' @keywords internal
#' @export
compute.ife2 <- function(formla, xformla, data, idname, tname, treat.t) {
    yname <- BMisc::lhs.vars(formla)
    dname <- BMisc::rhs.vars(formla)
    data$y <- data[,BMisc::lhs.vars(formla)]
    data$d <- data[,BMisc::rhs.vars(formla)]
    ##figure out the dates and make balanced panel
    tlist <- unique(data[,tname])[order(unique(data[,tname]))] ## this is going to be from smallest to largest

    ##flist <- unique(data[,first.treat.name])[order(unique(data[,first.treat.name]))]
    ##flist <- flist[flist>0]

    ##################################
    ## do some error checking
    if (!is.numeric(tlist)) {
        warning("not guaranteed to order time periods correclty if they are not numeric")
    }

    if (!(length(dname)==1)) {
        stop("you should only include one treatment variable")
    }

    ## check that all treatments are either 0 or 1
    ####################################

    
    tlen <- length(tlist)
}
