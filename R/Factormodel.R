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
    difft <- tlist[-c(1,2)] ## drop the first two time periods
    ## take X's from the first time period; if they do not vary over time this
    ## is fine...
    X <- model.matrix(xformla, data=data[ data[,tname]==tlist[1], ])
    first2pers <- BMisc::panel2cs(data[ data[,tname]==tlist[1] | data[,tname]==tlist[2], ], timevars="y", idname=idname, tname=tname)
    dY2 <- first2pers$dy
    XX <- cbind(X, dY2)
    Y1 <- data[ data[,tname]==tlist[1], ]$y

    Ytime <- do.call(cbind,lapply(tlist[-c(1,2)], function(tval) {
        data[ data[,tname]==tval,]$y
    })) ## matrix that contains outcomes in periods 3:(T-2)
    p <- mean(data$d)

    tt <- data[ data[,tname]==tlist[1], ]$d == 1## indicator for whether or not individual is in treated group (of length n)
    D <- 1*tt
    p <- mean(D)
    Xt <- (D/p)*X ##X[tt,]
    dY2t <- (D/p)*dY2 ##dY2[tt]
    Y1t <- (D/p)*Y1 ##Y1[tt]
    Ytimet <- (D/p)*Ytime ##Ytime[tt,]
    Xu <- ((1-D)/(1-p))*X ##X[(1-tt),]
    dY2u <- ((1-D)/(1-p))*dY2 ##dY2[(1-tt)]
    Y1u <- ((1-D)/(1-p))*Y1 ##Y1[(1-tt)]
    Ytimeu <- ((1-D)/(1-p))*Ytime ##Ytime[(1-tt),]
    nt <- sum(tt)
    nu <- sum((1-tt))
    n <- nrow(X)
    K <- ncol(X)
    
    ## set up the moment conditions to estimate with GMM
    ## params should be (T-2)x(k+1) where k is the dimension of X;
    ## it should have bet3 as first k elements, bet4 as next k, so on,
    ## last (T-2) elements should be F
    ## d will just be a dummy data.frame with n rows
    k <- ncol(X)
    g <- function(params, d) {
        bet <- matrix(params[1:((tlen-2)*k)], nrow=(tlen-2), byrow=TRUE) ## (T-2)xk matrix
        browser()

        ## see relevant sections in appendix of paper; this mimics notation there
        A0 <- kronecker(diag((tlen-2)), cbind(X,dY2))  #Xu,dY2u
        tstarlen <- length(which(tlist <= treat.t))
        A1a <- kronecker(diag((tstarlen-3)), cbind(X, dY2)) #Xt,dY2t
        ##A1b <- kronecker(diag((1+tlen-tstarlen)), matrix(data=0, nrow=n, ncol=(K+1)))
        A1 <- cbind(A1a, matrix(0, nrow=nrow(A1a), ncol=((1+tlen-tstarlen)*(K+1))))
        i1 <- D/p
        ii <- rep(1,length(i1))
        A2 <- kronecker(diag((tlen-1 + K)),ii)
        A <- magic::adiag(rbind(A0,A1),A2)

        ## now get the matrix of instruments in the serially uncorrelated case
        ZZinneru <- lapply(1:(tlen-2), function(i) {cbind(Xu,Ytimeu[,-i]) })
        ZZu <- do.call(magic::adiag, ZZinneru)
        ZZinnert <- lapply(1:(tstarlen-3), function(i) {cbind(Xt,Ytimet[,-i]) }) ## make sure this works in case with 4 periods
        ZZta <- do.call(magic::adiag, ZZinnert)
        ZZt <- cbind(ZZta, matrix(0, nrow=nrow(ZZta), ncol=((1+tlen-tstarlen)*(K+ncol(Ytimet)-1)))) ## add extra columns of 0s due to fewer available periods for treated group
        B2 <- kronecker(diag((tlen-1+K)),i1)
        B <- magic::adiag(rbind(ZZu,ZZt),2) ## for the directly estimated parameters A2 is the same for "regressors" and "instruments"

        ## vector of outcomes
        y <- c(Ytime-Y1,  ##already will select treated or untreated outcomes based on matrix B
               Ytime[,1:(tstarlen-3)]-Y1, ## outcomes for treated
               cbind(dY2, Ytime)-Y1, ## these are for average change in outcomes over time for treated group
               c(X))##do.call(magic::adiag,lapply(1:ncol(X), function(i) as.matrix(X[,i])))
        
        solve(t(B)%*%A%*%t(A)%*%B)%*%t(B)%*%A%*%solve(A)%*%y
        solve(t(A)%*%A)
        
        
        ## if we get singular matrix, use less efficient estimator
        if (abs(det(t(B)%*%A%*%t(A)%*%B)) < .Machine$double.eps) {
            warning("required matrix not invertible...\ntrying simpler (less efficient) approach...")
            A0 <- cbind(X, dY2) ## Xu, dY2u
            A2 <- kronecker(diag(2+K),ii)
            A <- magic::adiag(A0,A2)
            B0 <- cbind(Xu, Ytimeu[,1])
            B2 <- kronecker(diag(2+K),i1)
            B <- magic::adiag(B0,B2)
            delt <- matrix(nrow=((K+1)+2+K), ncol=(tlen-2))
            for (i in 1:(tlen-2)) {
                y <- c(Ytime[,i] - Y1, Ytime[,i]-Y1, dY2, c(X))
                delt[,i] <- solve(t(A)%*%B%*%t(B)%*%A)%*%t(A)%*%B%*%t(B)%*%y ## this is right,
                ## a good way to check is to confirm that moments of X are correect
            }
        }

        ## this works for getting beta_t and F_t, I think
        solve(t(XXu)%*%ZZu%*%t(ZZu)%*%XXu)%*%t(XXu)%*%ZZu%*%t(ZZu)%*%c(Ytimeu)
        
        

        ## but still can use moment conditions from treated group
        ## plus, should we estimate ATT in the same step
        
        
        ## get the "X" moments
        Ft <- tail(params, (tlen-2))
        mX1 <- t(Xt)%*%Xt%*%t(bet)/nt
        mX2 <- t(Xt)%*%(Ytimet-Y1t)/nt
        mX3 <- t(Xt)%*%outer(dY2t, Ft)/nt
        mX <- c(mX2-mX1-mX3)

        ## get the "Y" moments
        mY1 <- t(Ytimet)%*%Xt%*%t(bet)/nt
        mY2 <- t(Ytimet)%*%(Ytimet-Y1t)/nt
        mY3 <- t(Ytimet)%*%outer(dY2t, Ft)/nt
        ## here have to "throw away" moments where s=t
        mY <- mY2-mY1-mY3
        mY[upper.tri(mY) | lower.tri(mY)]
        mY <- matrix(mY[upper.tri(mY) | lower.tri(mY)], ncol=ncol(mY1))
        mY <- c(mY)

        ## still need to put in available moments from the untreated group,
        ## but this should get us going
        c(mX,mY)
    }

    n <- nrow(X)
    g(seq(1,(tlen-2)*(k+1)), data.frame(c=seq(1,n)))
    gmm(g=g, x=data.frame(c=seq(1,n)))
}
