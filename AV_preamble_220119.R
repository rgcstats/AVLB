library(sampling)
library(splines)
library(parallel)
library(biglm)

# utility function

rescale <- function(x,lower,upper){
  # creates a vector proportional to x truncated above and below
  #  which still sums to the same value as sum(x)
  # check and set search interval
  N <- length(x)
  if(any(x<=0)) stop("error: all values of x must be strictly positive")
  if((upper*N)<sum(x)) stop("even upper value sums to less than sum(x)")
  if((lower*N)>sum(x)) stop("even lower value sums to more than sum(x)")
  interval1 <- lower/max(x) # lambda is so small that even the biggest value gets truncated below
  interval2 <- upper/min(x) # lambda is so big that even the smallest value gets truncated above
  out <- pmax(pmin(x,upper),lower)
  optfn <- function(lambda){
    newx <- pmax(pmin(lambda*x,upper),lower)
    sum(newx)-sum(x)
  }
  lambda <- uniroot(optfn,interval=c(interval1,interval2))$root
  return(pmax(pmin(lambda*x,upper),lower))
}

scenarios1 <- expand.grid(n=c(500,1500),N=6000,R=5000,
                          beta=4,bw=c(0.5,1,2,5),
                          power.selnprob=c(0.5,1,2),maxratio.z=c(1,1.5,2.5),
                          true.var.index=1,work.var.index=c(1,2),
                          maxratio.selnprob=40,
                          mean.ty=NA,UB=NA)
scenarios2 <- expand.grid(n=1500,N=6000,R=5000,
                          beta=4,bw=c(0.5,1,2,5),
                          power.selnprob=0.5,maxratio.z=seq(from=1,to=4,by=0.25),
                          true.var.index=1,work.var.index=1,
                          maxratio.selnprob=40,
                          mean.ty=NA,UB=NA)
scenarios3 <- expand.grid(n=25e3,N=100e3,R=15000,
                          beta=4,bw=c(0.5,1,2,5),
                          power.selnprob=c(0.5,1,2),maxratio.z=c(1,1.5,2.5),
                          true.var.index=1,work.var.index=1,
                          maxratio.selnprob=40,
                          mean.ty=NA,UB=NA)
scenarios <- rbind(scenarios1,scenarios2,scenarios3)
scenarios <- scenarios[order(scenarios$n,scenarios$bw,scenarios$power.selnprob,
                             scenarios$maxratio.z,scenarios$work.var.index),]
dups <- NULL
for(k in c(2:nrow(scenarios))){
  if(all(scenarios[k,]==scenarios[k-1,],na.rm=TRUE)) dups <- c(dups,k)
}
scenarios <- scenarios[-dups,]
# only use work.var.index=2 when n=500
scenarios <- scenarios[!((scenarios$n==1500)&(scenarios$work.var.index==2)),]

K <- nrow(scenarios)
scenarios$k <- 1:K

simulate.fn <- function(k,Rdiv,scenarios){
  R <- scenarios$R[k]/Rdiv
  simresults <- data.frame(k=rep(k,R*100),r=NA,method="none",est=NA,AIC=NA,true=NA,
                           model.chosen=NA,stringsAsFactors=FALSE)
  i <- 1
  set.seed(51850175)
  cat("Starting Scenario ",k," \n")
  # load parameters for this scenario
  n <- scenarios$n[k]
  N <- scenarios$N[k]
  beta <- scenarios$beta[k]
  maxratio.selnprob <- scenarios$maxratio.selnprob[k]
  power.selnprob <- scenarios$power.selnprob[k]
  bw <- scenarios$bw[k]
  maxratio.z <- scenarios$maxratio.z[k]
  true.var.index <- scenarios$true.var.index[k]
  work.var.index <- scenarios$work.var.index[k]
  # generate population
  #xvals <- c(1:100)/100
  xvals <- exp(qnorm(p=c(1:100)/101,mean=-1/32,sd=0.25))
  xvals.group <- rep(c(1:5),each=20)
  popn <- data.frame(x=rep(xvals,each=N/100),xgroup=rep(xvals.group,each=N/100),
                     z=rep(c(0,1),times=N/2))
  popn$xgrp.mean <- ave(popn$x,popn$xgroup,FUN=mean)
  # make probs of selection proportional to sqrt(x)? triple for industry 1
  # put in  4 steps into the probabilities of selection, make it basically on
  #   sqrt(x), jump up by 0.1 at each cutpoint
  popn$selnprob <- popn$x^power.selnprob * (1+(maxratio.z-1)*popn$z)
  popn$selnprob <- popn$selnprob / sum(popn$selnprob) * n
  popn$selnprob <- rescale(popn$selnprob,1/maxratio.selnprob,1)
  # calculate UB
  #popn$var.y.given.x <- popn$mu^2/shape.y
  popn$var.y.given.x <- 0.25*popn$x^true.var.index
  scenarios$UB[k] <- sum((1/popn$selnprob-1)*popn$var.y.given.x)
  allcoef <- NULL
  # now simulate R y-popns and samples with predictions of ty
  for(r in c(1:R)){
    simresults$r[i:(nrow(simresults))] <- r
    # generate population y values
    popn$mu <- beta*popn$x + sin(popn$x/bw*2*pi)
    popn$y <- popn$mu + rnorm(N)*sqrt(popn$var.y.given.x)
    simresults[i:nrow(simresults),"true"] <- sum(popn$y)
    popn <- popn[rank(runif(N)),]
    s <- UPsystematic(popn$selnprob)
    samp <- popn[s==1,]
    nonsamp <- popn[s==0,]
    # truncate nonsamp and popn x values to avoid extrapolation problems
    trunc.nonsamp.x <- pmax(pmin(nonsamp$x,max(samp$x)),min(samp$x))
    trunc.popn.x <- pmax(pmin(popn$x,max(samp$x)),min(samp$x))
    # HT estimator
    simresults[i,"method"] <- "HT"
    simresults[i,"est"] <- sum(samp$y/samp$selnprob)
    i <- i+1
    # fit lm model
    lmfit0 <- lm(y~1,data=samp,weights=1/samp$x^work.var.index)
    lmfitz0 <- lm(y~z,data=samp,weights=1/samp$x^work.var.index)
    lmfit <- lm(y~x,data=samp,weights=1/samp$x^work.var.index)
    allcoef <- rbind(allcoef,lmfit$coefficients)
    lmfitz <- lm(y~x+z,data=samp,weights=1/samp$x^work.var.index)
    # select spline or lm models based on AIC
    # models: 1=intercept, 2=z, 3=x, 4=x+z
    #         5=no knots, 6=no knots+z,
    #         7,9,11,13,15,17,19,21,23,25 = spline with 1-10 knots,
    #         8,10,12,14,16,18,20,22,24,26 = z + spline with 1-10 knots
    model.chosen <- 1 ; bestmodel <- lmfit0
    if(BIC(lmfitz0)<BIC(bestmodel)){
      model.chosen <- 2 ; bestmodel <- lmfitz0
    }
    if(BIC(lmfit)<BIC(bestmodel)){
      model.chosen <- 3 ; bestmodel <- lmfit
    }
    if(BIC(lmfitz)<BIC(bestmodel)){
      model.chosen <- 4 ; bestmodel <- lmfitz
    }
    candidate <- lm(y~bs(x,knots=NULL,Boundary.knots=c(0,max(popn$x)),degree=3),
                    data=samp,weights=1/samp$x^work.var.index)
    if(BIC(candidate)<BIC(bestmodel)){
      model.chosen <- 5 ; bestmodel <- candidate
    }
    candidate <- lm(y~z+bs(x,knots=NULL,Boundary.knots=c(0,max(popn$x)),degree=3),
                         data=samp,weights=1/samp$x^work.var.index)
    if(BIC(candidate)<BIC(bestmodel)){
      model.chosen <- 6 ; bestmodel <- candidate
    }
    j <- 7
    for(numknots in c(1:10)){
      candidate <- lm(y~bs(x,knots=min(popn$x)+c(1:numknots)/(numknots+1)*(max(popn$x)-min(popn$x)),Boundary.knots=c(min(popn$x),max(popn$x)),degree=3),data=samp,weights=1/samp$x^work.var.index)
      if(BIC(candidate)<BIC(bestmodel)){
        model.chosen <- j ; bestmodel <- candidate
      }
      j <- j+1
      candidate <- lm(y~z+bs(x,knots=min(popn$x)+c(1:numknots)/(numknots+1)*(max(popn$x)-min(popn$x)),Boundary.knots=c(min(popn$x),max(popn$x)),degree=3),data=samp,weights=1/samp$x^work.var.index)
      if(BIC(candidate)<BIC(bestmodel)){
        model.chosen <- j ; bestmodel <- candidate
      }
      j <- j+1
    }
    # calculate ests of ty from linear model fits
    simresults[i,"method"] <- "lm"
    simresults[i,"est"] <- sum(samp$y)+sum(predict(lmfit,newdata=nonsamp))
    i <- i+1
    # calculate GREG and BLUP from best model
    simresults[i,"method"] <- "best"
    simresults[i,"est"] <- sum(samp$y)+sum(predict(bestmodel,newdata=data.frame(x=trunc.nonsamp.x,z=nonsamp$z)))
    simresults[i,"model.chosen"] <- model.chosen
    i <- i+1
  }
  simresults <- simresults[simresults$method!="none",]
  simresults <- merge(simresults,scenarios,by="k")
  simresults$error <- simresults$est - simresults$true
  simresults
}

