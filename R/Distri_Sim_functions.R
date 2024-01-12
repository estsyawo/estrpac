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





#===============================================================================================>
#' \code{BB_CBOM} computes Montiel Olea & Plagborg-Moller 2018 Bayesian simultaneous bands.
#'
#' @param DQmat G by L matrix of random draws from the joint distribution
#' @param alpha nominal level for \eqn{100(1-\alpha)\%} simultaneous coverage
#' @param method type of credible interval to calibrate for simultaneous coverage. Defaults to
#' the equal-tailed interval (ETI); other options include  'HDI', 'BCI' and 'SI'. 
#' See \link[bayestestR]{ci} for details.
#' @param qrangeID an index of the rows of \code{DQmat} to obtain simultaneous coverage for.
#' @param tol a tolerance level as stopping rule for the algorithm.
#' 
#' @return BB: Lower and upper Bayesian bands
#' @return zeta: the calibrated level at which simultaneous coverage is achieved.
#' 
#' @importFrom stats quantile
#' @importFrom utils tail
#' @importFrom bayestestR ci
#' 
#' @examples 
#' set.seed(1); BB_CBOM(matrix(rnorm(1000),nrow=5))
#' @export

BB_CBOM=function (DQmat, alpha = 0.1,method = "ETI",qrangeID=NULL,tol=.Machine$double.eps^0.25) 
{
  G = nrow(DQmat); L = ncol(DQmat)
  if(is.null(qrangeID)){qrangeID=1:G}
  zbnds=c((alpha/G), alpha)/2
  # function to compute credible intervals
  cif=function(g,alpha) {
    if(method!="ETI"){
      ciobj=ci(DQmat[g,],ci=(1-alpha),method=method)
      ans=c(ciobj$CI_low,ciobj$CI_high)
    }else{ans = c(quantile(DQmat[g,],probs = c(alpha/2,(1-alpha/2)), type = 1))}
    ans
  }
  
  fn = function(zeta){
    vm = t(sapply(1:G,cif,alpha=zeta*2))
    zf = function(l){all(DQmat[qrangeID,l]>=vm[qrangeID,1]) & 
        all(DQmat[qrangeID,l]<=vm[qrangeID,2])}
    as.numeric(sum(sapply(1:L, zf))>=(1-alpha)*L)
  } 
  fn=Vectorize(fn); zb = zbnds; sq = seq(zb[1],zb[2],length.out = 10)
  fnsq = fn(sq); dl = sq[2]-sq[1]
  if(fnsq[1]==0){
    warning("Coverage may be conservative over the entire search interval")
    zeta=sq[1]
    BB = t(sapply(1:G,cif,alpha=zeta*2))
  }else if(tail(fnsq,1)==0 && fnsq[1]==1){ #search over interval
    while(dl>tol){#need dl<= tol
      d = which(fnsq==0)[1]
      zb=sq[c((d-1),d)]
      sq = seq(zb[1],zb[2],length.out = 10)
      fnsq=fn(sq)
      dl = sq[2]-sq[1]
      #print(dl)
    } #end while()
    #compute zeta and Bayesian Bands
    zeta=sq[which(fnsq==0)[1]-1]
    BB = t(sapply(1:G,cif,alpha=zeta*2))
  }else{
    warning("Coverage is 1-alpha over the entire search interval")
    zeta=tail(sq,1)
    BB = t(sapply(1:G,cif,alpha=zeta*2))
  }
  list(BB = BB, zeta = zeta)
}
#===============================================================================================>





#=====================================================================================>
#' \code{Band.Fn.objs} constructs Chernozhukov et al. 2013 Uniform Confidence Bands
#'
#' @param Fn G-length vector of the function of interest
#' @param Boot.Fn a \eqn{G\times L} matrix of samples of the function Fn
#' @param n sample size of the data used to compute Fn
#' 
#' @return Sigv:
#' @return tbv: 
#' 
#' @importFrom stats IQR qnorm
#' 
#' @examples 
#' set.seed(1); Band.Fn.objs(Fn=rep(0,5),Boot.Fn=matrix(rnorm(1000),nrow=5),n=100)
#' @export

Band.Fn.objs<- function(Fn,Boot.Fn,n){
  Zbstar = Boot.Fn
  for (j in 1:ncol(Boot.Fn)){Zbstar[,j] = Zbstar[,j] - Fn}
  Zbstar = sqrt(n)*Zbstar
  Sigv = apply(Zbstar,1,IQR)/(qnorm(0.75)-qnorm(0.25)) #a robust estimator of sigma
  
  tbmat = abs(Zbstar)
  for (j in 1:ncol(tbmat)) {tbmat[,j]/Sigv}
  
  tbv = apply(tbmat,2,max)
  list(Sigv=Sigv,tbv=tbv)
}
#=====================================================================================>
#' \code{UBand.Fn} constructs Chernozhukov et al. 2013 Uniform Confidence Bands
#' using output from \link{Band.Fn.objs}
#' @param Fn G-length vector of the function of interest
#' @param BDobjs output object from \link{Band.Fn.objs}
#' @param n sample size of the data used to compute Fn
#' @param alpha significance level
#' 
#' @return LB.Fn: G-length lower Uniform Confidence band
#' @return UB.Fn: G-length upper Uniform Confidence band
#' 
#' @importFrom stats quantile
#' 
#' @examples 
#' set.seed(1); 
#' BDobjs=Band.Fn.objs(Fn=rep(0,5),Boot.Fn=matrix(rnorm(1000),nrow=5),n=100)
#' UBand.Fn(Fn=rep(0,5),BDobjs=BDobjs,n=100)
#' 
#' @export

UBand.Fn=function(Fn,BDobjs,n,alpha=0.05){
  LB.Fn = Fn - quantile(BDobjs$tbv,1-alpha)*BDobjs$Sigv/sqrt(n) 
  UB.Fn = Fn + quantile(BDobjs$tbv,1-alpha)*BDobjs$Sigv/sqrt(n)
  list(LB.Fn=LB.Fn, UB.Fn=UB.Fn)
}
#=====================================================================================>
#' \code{compute.Ubands} is a high-level function for computing uniform confidence bands 
#' according to Chernozhukov et al. 2013
#' @param Fn G-length vector of the function of interest
#' @param Boot.Fn a \eqn{G \times L} matrix of samples of the function Fn
#' @param n sample size of the data used to compute Fn
#' @param alpha significance level
#' 
#' @return LB.Fn: G-length lower Uniform Confidence band
#' @return UB.Fn: G-length upper Uniform Confidence band
#' 
#' @importFrom stats quantile
#' 
#' @examples 
#' set.seed(1); 
#' Ubands=compute.Ubands(Fn=rep(0,5),Boot.Fn=matrix(rnorm(1000),nrow=5),n=100)
#' 
#' @export

compute.Ubands<- function(Fn,Boot.Fn,n,alpha = 0.05){
  BandObjects = Band.Fn.objs(Fn=Fn,Boot.Fn = Boot.Fn,n=n)
  UBand.Fn(Fn=Fn,BDobjs = BandObjects,alpha = alpha,n=n)
}

#=====================================================================================>
