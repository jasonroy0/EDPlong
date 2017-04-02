## here we place functions written in R

expit <- function(x) 1/(1+exp(x))
inv<-function(x) chol2inv(chol(x))

dinvchisq <- function (
  ### density function and for the (scaled) inverse-chi-squared distribution.
  x    			##<< vector of quantiles
  , df				##<< degrees of freedom parameter, usually represented as nu
  , scale = 1/df		##<< scale parameter, usually represented as lambda.
  , log = FALSE		##<< Logical. If log=TRUE, then the logarithm of the density is returned.
){
  # adopted from copyLeft pakcage LaplaceDemon, as not available for R 2.10 which runs on the cluster
  ##seealso<< \code{\link{dinvchisq}}
  x <- as.vector(x)
  df <- as.vector(df)
  scale <- as.vector(scale)
  if (any(x <= 0)) 
    stop("x must be positive.")
  if (any(df <= 0)) 
    stop("The df parameter must be positive.")
  if (any(scale <= 0)) 
    stop("The scale parameter must be positive.")
  NN <- max(length(x), length(df), length(scale))
  x <- rep(x, len = NN)
  df <- rep(df, len = NN)
  scale <- rep(scale, len = NN)
  nu <- df/2
  dens <- nu * log(nu) - log(gamma(nu)) + nu * log(scale) - 
    (nu + 1) * log(x) - (nu * scale/x)
  if (log == FALSE) 
    dens <- exp(dens)
  dens
}
attr(dinvchisq,"ex") <- function(){
  x <- dinvchisq(1,1,1)
  x <- rinvchisq(10,1)
  
  #Plot Probability Functions
  x <- seq(from=0.1, to=5, by=0.01)
  plot(x, dinvchisq(x,0.5,1), ylim=c(0,1), type="l", main="Probability Function",
       ylab="density", col="red")
  lines(x, dinvchisq(x,1,1), type="l", col="green")
  lines(x, dinvchisq(x,5,1), type="l", col="blue")
  legend(3, 0.9, expression(paste(nu==0.5, ", ", lambda==1),
                            paste(nu==1, ", ", lambda==1), paste(nu==5, ", ", lambda==1)),
         lty=c(1,1,1), col=c("red","green","blue"))	
}


rinvchisq <- function (
  ### random number function and for the (scaled) inverse-chi-squared distribution.
  n				## the number of observations. If length(n) > 1, then the length is taken to be the number required.
  , df
  , scale = 1/df
){
  ##seealso<< \code{\link{rinvchisq}}
  df <- rep(df, len = n)
  scale <- rep(scale, len = n)
  if (any(df <= 0)) 
    stop("The df parameter must be positive.")
  if (any(scale <= 0)) 
    stop("The scale parameter must be positive.")
  x <- (df * scale)/rchisq(n, df = df)
  return(x)
}

############## REGRESSION UPDATES ######################################
#regression variance update
newregvar<-function(tempx,tempy){
  an<-beta.a0+length(tempy)/2
  ifelse((ncol(as.matrix(tempx))>1),(xtx<-t(tempx)%*%tempx),(xtx<-as.matrix(tempx)%*%t(as.matrix(tempx))))
  newprec<-xtx+prec0
  if(ncol(as.matrix(tempx))>1){betan<-(inv(xtx+prec0))%*%(prec0%*%beta0+t(tempx)%*%tempy)}
  if(ncol(as.matrix(tempx))<=1){betan<-(inv(xtx+prec0))%*%(prec0%*%beta0+(tempx*tempy))}
  bn<-beta.b0+.5*(crossprod(tempy,tempy)+t(beta0)%*%prec0%*%beta0-t(betan)%*%newprec%*%betan)
  return(rinvgamma(1,an,bn))
}

#regression betas update
newbet<-function(tempx,tempy,sig){
  ifelse((ncol(as.matrix(tempx))>1),(xtx<-t(tempx)%*%tempx),(xtx<-as.matrix(tempx)%*%t(as.matrix(tempx))))
  newprec<-xtx+prec0
  if(ncol(as.matrix(tempx))>1){betan<-(inv(xtx+prec0))%*%(prec0%*%beta0+t(tempx)%*%tempy)}
  if(ncol(as.matrix(tempx))<=1){betan<-(inv(xtx+prec0))%*%(prec0%*%beta0+(tempx*tempy))}
  return(rmvn(1,betan,sig*inv(newprec)))
}
###########################################################################


samp.b<-function(tempz,tempy,sigy,sigb){
  ifelse((ncol(as.matrix(tempz))>1),(ztz <- t(tempz)%*%tempz) , ztz <- outer(tempz,tempz)) 
  new.sig <- inv(ztz/sigy+diag(nZ)/sigb)
  ifelse((ncol(as.matrix(tempz))>1),new.mu <- new.sig%*%t(tempz)%*%tempy / sigy , new.mu <- new.sig%*%tempz%*%tempy / sigy) 
  return(rmvn(1,new.mu,new.sig))
}

samp.u<- function(tempy,sigy,sigu) {
  tempn <- length(tempy)
  new.sig <- (tempn/sigy + 1/sigu)^-1
  new.mu <- new.sig * sum(tempy)/sigy
  return(rnorm(1,new.mu,sqrt(new.sig)))
}

## see https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf (section 5) for more details

############CONTINUOUS COVARIATES UPDATES##################################
updatevar<-function(tempx){
  newdf<-nu0+length(tempx)
  if(length(tempx)==1){varx<-0}
  if(length(tempx)>1){varx<-var(tempx)}
  numer<-nu0*tau0+(length(tempx)-1)*varx+(c0*length(tempx)/(c0+length(tempx)))*(mean(tempx)-mu0)^2
  newval<-rinvchisq(1,newdf,numer/newdf)
  return(newval)
}

updatemean<-function(tempx,curtau){
  newvar<-1/(c0/curtau+length(tempx)/curtau)
  newmean<-(mu0*c0/curtau+mean(tempx)*length(tempx)/curtau)*newvar
  newval<-rnorm(1,newmean,sqrt(newvar))
  return(newval)
}
###########################################################################



##############MASS PARAMETER UPDATES#######################################
newalpthet<-function(numclus,prioralp,alpa0,alpb0,nobs){
  eta<-rbeta(1,prioralp+1,nobs)
  pieta<-(numclus/(nobs*(1-log(eta))))/(1+numclus/(nobs*(1-log(eta))))
  whichmix<-rbinom(1,1,pieta)
  newalp<-whichmix*rgamma(1,(alpa0+numclus),(alpb0-log(eta)))+
    (1-whichmix)*rgamma(1,(alpa0+numclus-1),(alpb0-log(eta)))
  return(newalp)
}

#update alpha_psi using MH
newalppsi<-function(alpha,cluster,alpa0,alpb0){
  sortuniquey<-sort(unique(cluster[,1]))
  ss<-numeric(length(sortuniquey))
  for(j in 1:length(sortuniquey)){
    ss[j]<-sum(cluster[,1]==sortuniquey[j])
  }
  likecur<-dgamma(alpha,alpa0,alpb0)*(alpha^(nrow(unique(s))-k))*prod((alpha+ss)*beta(alpha+1,ss))
  #proposed<-rgamma(1,2,1)
  proposed <- rgamma(1,1,3)
  likeprop<-dgamma(proposed,alpa0,alpb0)*(proposed^(nrow(unique(s))-k))*prod((proposed+ss)*beta(proposed+1,ss))
  rat<-(likeprop)/(likecur)
  if(runif(1,0,1)<rat){alpha<-proposed}
  return(alpha)
}
############################################################################
