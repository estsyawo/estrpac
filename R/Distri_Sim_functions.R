#===========================================================================================#
#' A Generic Independence Metropolis-Hastings Algorithm
#'
#' \code{IndepMHgen} computes random draws of parameters using a normal proposal distribution.
#' This function implements a generic form of \code{IndepMH} from the package \code{bayesdistreg}
#'
#' @param start starting values of parameters for the MH algorithm.
#' It is automatically generated from the normal proposal distribution but the user can also specify.
#' @param posterior the log posterior distribution function.
#' Should take parameter input of the same length as \code{start} or \code{propob$mode} 
#' @param ... additional arguments to the posterior function
#' @param propob a list of mode and variance-covariance matrix of the normal proposal distribution. 
#' Save list as propob=list(mode=mode,var=variance-covariance)
#' @param const a vector function of parameters showing non-negative inequality constraints to be satisfied. 
#' @param seed an integer as seed for reproducibility
#' @param scale a value multiplied by \code{propob$var} in order to adjust the proposal distribution. Else
#' set to the character string "HS18" for the Herbst and Shorfheide (2018) scale updating for a 0.25
#' acceptance ratio.
#' The default is \code{1.5} but the user can adjust it until a satisfactory acceptance rate is obtained.
#' @param iter number of random draws desired (default: 5000)
#' @param burn burn-in period for the MH algorithm (default: floor(iter/10))
#' @param report a numeric frequency (i.e. after how many iterations to report progress) for reporting
#'  algorithm progress; default - NULL
#' @return Matpram a matrix of parameter draws
#' @return postvals vector of posterior values corresponding to parameter draws \code{Matpram}
#' @return AcceptRatio the acceptance ratio
#' 
#' @examples 
#' #a toy example for illustration
#' ## f(c) = 1/(3.618*sqrt(pi))* exp(-0.6*(c[1]-2)^2-0.4*(c[2]+2)^2) 
#' # an improper posterior
#' logpost = function(c) -0.6*(c[1]-2)^2-0.4*(c[2]+2)^2 #log posterior distribution
#' optp<-optim(par=c(0,0),fn=logpost,control=list(fnscale=-1),hessian = TRUE) 
#' # laplace approximation of the posterior
#' propob = list(mode=optp$par,var=-solve(optp$hessian)) #parameters of proposal distribution
#' eigen(propob$var)$values # var-cov of proposal distribution is positive definite
#' MHobj<- indepMHgen(posterior = logpost,propob = propob,scale = "HS18",iter = 6000,report=30)
#' # create an independent Metropolis-Hastings object
#' dim(MHobj$Matpram) # a 2 x 5000 matrix with columns corresponding to draws of c1 and c2
#' par(mfrow=c(1,2))
#' hist(MHobj$Matpram[1,],20,main = "Histogram c1",xlab = "c1")
#' hist(MHobj$Matpram[2,],20,main = "Histogram c2",xlab = "c2"); par(mfrow=c(1,1))
#' MHobj$AcceptRatio # acceptance ratio
#'
#' @export

indepMHgen<- function(start=NULL,posterior=NULL,...,propob=NULL,const=NULL,
                      seed=1,scale=1.5,iter=5000,burn=floor(0.1*iter),
                      report=NULL){
  if(scale!="HS18"){
    varprop = scale*propob$var
  }else{
    varprop = propob$var
    c0 = 1.0 #initialise adapting scale
    }
  if(!is.null(seed)){set.seed(seed = seed)}
  npar = length(propob$mode)
  Mat = array(0, c(iter, npar)); postvals<- c(0)

  if(is.null(const))
  {  
    if(is.null(start)){
      start = MASS::mvrnorm(n=1,propob$mode,varprop)
    }
    Mat[1,] = start; AccptRate<-0; postvals[1]<- posterior(start,...)
    
    for(i in 2:iter){
    start= Mat[i-1,]
    prop = MASS::mvrnorm(n=1,propob$mode,varprop)#make a draw from proposal dist
    lpa = posterior(prop,...); lpb = postvals[i-1]
    accprob = exp(lpa-lpb)
    if(is.na(accprob)){accprob=0} #penalise NA fun values
    # the other part cancels out because the normal distribution is symmetric
    if(stats::runif(1)< accprob){
      Mat[i,]=prop
      AccptRate<- AccptRate +1
      postvals[i]<- lpa
    }else{
      Mat[i,]=start 
      postvals[i]<- postvals[i-1]
    }
    if(scale=="HS18"){
      rx = AccptRate/i #acceptance ratio so far
      c1 = c0*(0.95 + 0.1*exp(16*(rx-0.25))/(1+exp(16*(rx-0.25))))
      c0 = c1
    }
    if(!is.null(report)){
      if(i%%report==0){message(i," iterations done. ",I(iter-i)," more to go. \n")}
    }
  }
  }else{
    if(is.null(start)){
      cnt<- 0
      start = MASS::mvrnorm(n=1,propob$mode,varprop)
      while(any(const(start)<0)){
        start = MASS::mvrnorm(n=1,propob$mode,varprop) ; cnt<- cnt+1
        if(cnt>burn){stop("Cannot generate proposal draws within the constraint region")}
      }
    }
    Mat[1,] = start; AccptRate<-0; postvals[1]<- posterior(start,...)
    
    for(i in 2:iter){
      start= Mat[i-1,]
      prop = MASS::mvrnorm(n=1,propob$mode,varprop)#make a draw from proposal dist
      cnt<- 0
      while(any(const(prop)<0)){
        prop = MASS::mvrnorm(n=1,propob$mode,varprop) ; cnt<- cnt+1
        if(cnt>burn){stop("Cannot generate proposal draws within the constraint region")}
      }
      lpa = posterior(prop,...); lpb = postvals[i-1]
      accprob = exp(lpa-lpb)
      # the other part cancels out because the normal distribution is symmetric
      if(is.na(accprob)){accprob=0} #penalise NA fun values
      if(stats::runif(1)< accprob){
        Mat[i,]=prop
        AccptRate<- AccptRate +1
        postvals[i]<- lpa
      }else{
        Mat[i,]=start 
        postvals[i]<- postvals[i-1]
      }
      if(scale=="HS18"){
        rx = AccptRate/i #acceptance ratio so far
        c1 = c0*(0.95 + 0.1*exp(16*(rx-0.25))/(1+exp(16*(rx-0.25))))
        c0 = c1
      }
      if(!is.null(report)){
        if(i%%report==0){message(i," iterations done. ",I(iter-i)," more to go. \n")}
      }
    }
  }
  if(!is.null(report)){message("indepMHgen algorithm successful\n")}
  val = list(Matpram=t(Mat[-c(1:burn),]),postvals=postvals[-c(1:burn)],AcceptRatio = AccptRate/iter)
  return(val)
}
#===========================================================================================#


#===========================================================================================#
#' A Generic Random Walk Metropolis-Hastings Algorithm
#'
#' \code{rwMHgen} computes random draws of parameters using a normal proposal distribution.
#' This function implements a generic form of \code{RWMH} from the package \code{bayesdistreg}
#'
#' @param start starting values of parameters for the MH algorithm.
#' It is automatically generated from the normal proposal distribution but the user can also specify.
#' @param posterior the log posterior distribution function. 
#' Should take parameter input of the same length as \code{start} or \code{propob$mode}
#' @param ... additional arguments to the posterior function
#' @param propob a list of mode and variance-covariance matrix of the normal proposal distribution. 
#' Save list as propob=list(mode=mode,var=variance-covariance)
#' @param const a vector function of parameters showing non-negative inequality constraints to be satisfied. 
#' @param seed an integer as seed for reproducibility
#' @param scale a value multiplied by \code{propob$var} in order to adjust the proposal distribution. Else
#' set to the character string "HS18" for the Herbst and Shorfheide (2018) scale updating for a 0.25
#' acceptance ratio.
#' The default is \code{1.5} but the user can adjust it until a satisfactory acceptance rate is obtained.
#' @param iter number of random draws desired (default: 5000)
#' @param burn burn-in period for the MH algorithm (default: floor(0.1*iter))
#' @param report a numeric frequency (i.e. after how many iterations to report progress) for reporting
#'  algorithm progress; default - NULL
#' @return Matpram a matrix of parameter draws
#' @return postvals vector of posterior values corresponding to parameter draws \code{Matpram}
#' @return AcceptRatio the acceptance ratio
#' 
#' @examples 
#' #a toy example for illustration
#' ## f(c) = 1/(3.618*sqrt(pi))* exp(-0.6*(c[1]-2)^2-0.4*(c[2]+2)^2) 
#' # an improper posterior
#' logpost = function(c) -0.6*(c[1]-2)^2-0.4*(c[2]+2)^2 #log posterior distribution
#' optp<-optim(par=c(0,0),fn=logpost,control=list(fnscale=-1),hessian = TRUE) 
#' # laplace approximation of the posterior
#' propob = list(mode=optp$par,var=-solve(optp$hessian)) #parameters of proposal distribution
#' eigen(propob$var)$values # var-cov of proposal distribution is positive definite
#' MHobj<- rwMHgen(posterior = logpost,propob = propob,scale = "HS18",iter = 6000,report=20)
#' # create an independent Metropolis-Hastings object
#' dim(MHobj$Matpram) # a 2 x 5000 matrix with columns corresponding to draws of c1 and c2
#' par(mfrow=c(1,2))
#' hist(MHobj$Matpram[1,],20,main = "Histogram c1",xlab = "c1")
#' hist(MHobj$Matpram[2,],20,main = "Histogram c2",xlab = "c2"); par(mfrow=c(1,2))
#' MHobj$AcceptRatio # acceptance ratio
#'
#' @export

rwMHgen<- function(start=NULL,posterior=NULL,...,propob=NULL,const=NULL,
                      seed=1,scale=1.5,iter=5000,burn=floor(0.1*iter),
                      report=NULL){
  if(scale!="HS18"){
    varprop = scale*propob$var
  }else{
    varprop = propob$var
    c0 = 1.0 #initialise adapting scale
  }
  if(!is.null(seed)){set.seed(seed = seed)}
  npar = length(propob$mode)
  Mat = array(0, c(iter, npar)); postvals<- c(0)
  
  if(is.null(const))
  {  
    if(is.null(start)){
      start = MASS::mvrnorm(n=1,propob$mode,varprop)
    }
    Mat[1,] = start; AccptRate<-0; postvals[1]<- posterior(start,...)
    
    for(i in 2:iter){
      start= Mat[i-1,]
      prop = MASS::mvrnorm(n=1,start,varprop)#make a draw from proposal dist
      lpa = posterior(prop,...); lpb = postvals[i-1]
      accprob = exp(lpa-lpb)
      if(is.na(accprob)){accprob=0} #penalise NA fun values
      # the other part cancels out because the normal distribution is symmetric
      if(stats::runif(1)< accprob){
        Mat[i,]=prop
        AccptRate<- AccptRate +1
        postvals[i]<- lpa
      }else{
        Mat[i,]=start 
        postvals[i]<- postvals[i-1]
      }
      if(scale=="HS18"){
        rx = AccptRate/i #acceptance ratio so far
        c1 = c0*(0.95 + 0.1*exp(16*(rx-0.25))/(1+exp(16*(rx-0.25))))
        c0 = c1
      }
      if(!is.null(report)){
        if(i%%report==0){message(i," iterations done. ",I(iter-i)," more to go. \n")}
      }
    }
  }else{
    if(is.null(start)){
      cnt<- 0
      start = MASS::mvrnorm(n=1,propob$mode,varprop)
      while(any(const(start)<0)){
        start = MASS::mvrnorm(n=1,propob$mode,varprop) ; cnt<- cnt+1
        if(cnt>burn){stop("Cannot generate proposal draws within the constraint region")}
      }
    }
    Mat[1,] = start; AccptRate<-0; postvals[1]<- posterior(start,...)
    
    for(i in 2:iter){
      start= Mat[i-1,]
      prop = MASS::mvrnorm(n=1,start,varprop)#make a draw from proposal dist
      cnt<- 0
      while(any(const(prop)<0)){
        prop = MASS::mvrnorm(n=1,start,varprop) ; cnt<- cnt+1
        if(cnt>burn){stop("Cannot generate proposal draws within the constraint region")}
      }
      lpa = posterior(prop,...); lpb = postvals[i-1]
      accprob = exp(lpa-lpb)
      # the other part cancels out because the normal distribution is symmetric
      if(is.na(accprob)){accprob=0} #penalise NA fun values
      if(stats::runif(1)< accprob){
        Mat[i,]=prop
        AccptRate<- AccptRate +1
        postvals[i]<- lpa
      }else{
        Mat[i,]=start 
        postvals[i]<- postvals[i-1]
      }
      if(scale=="HS18"){
        rx = AccptRate/i #acceptance ratio so far
        c1 = c0*(0.95 + 0.1*exp(16*(rx-0.25))/(1+exp(16*(rx-0.25))))
        c0 = c1
      }
      if(!is.null(report)){
        if(i%%report==0){message(i," iterations done. ",I(iter-i)," more to go. \n")}
      }
    }
  }
  if(!is.null(report)){message("rwMHgen algorithm successful\n")}
  val = list(Matpram=t(Mat[-c(1:burn),]),postvals=postvals[-c(1:burn)],AcceptRatio = AccptRate/iter)
  return(val)
}
#===========================================================================================#