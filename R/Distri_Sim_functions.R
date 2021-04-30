#' @export
aclass<- function(dat,ndrs=NULL,seed=NULL,W=NULL){# a generic class creator
  if(is.null(W)){
    val = structure(list(dat=dat,ndrs=ndrs,seed=seed), class="default")
    return(val)#non-weighted bayesian normal reg as a class
  }else{
    val =  structure(list(dat=dat,ndrs=ndrs,seed=seed,W=W), class="weight")
    return(val)
  }
}



#===========================================================================================#
#' @export
simsigbet <- function(xx) UseMethod("simsigbet")

#===========================================================================================#

#' @importFrom invgamma rinvgamma
#' 
#' @export
simsigbet.default<- function(xx){
  dat=as.matrix(xx$dat); ndrs=xx$ndrs;seed=xx$seed
  if(is.null(ndrs)){
    ndrs=10000
  }
  if(is.null(seed)){
    seed=1
  }
  y = dat[,1]; x = as.matrix(dat[,-1]); nr=nrow(x)
  reg<-lm(y~as.matrix(x))
  bet<-reg$coefficients; nb=length(bet)
  Vb<- vcov(reg)/sigma(reg);
  S = sum(reg$residuals^2)
  set.seed(seed = seed)
  sig<- 1/stats::rgamma(ndrs,shape = (nr-nb)/2,scale = 2/S)
  obtbets<- function(sig){
    set.seed(seed)
    beta=MASS::mvrnorm(n=1,bet,sig*Vb)
    return(beta)
  }
  beta=sapply(sig,obtbets)

  val<- list(beta=t(beta),sig=sig)
  return(val)
}



#===========================================================================================#

#' @importFrom invgamma rinvgamma
#'
#' @export

simsigbet.weight<- function(xx){
  dat=as.matrix(xx$dat); ndrs=xx$ndrs;seed=xx$seed;W=xx$W
  if(is.null(ndrs)){
    ndrs=10000
  }
  if(is.null(seed)){
    seed=1
  }
  if(mean(W)<1){
    W = W/mean(W) # scale up to make mean W=1
  }
  y = dat[,1]; x = as.matrix(dat[,-1]); 
  reg<-lm(y~as.matrix(x),weights = W)
  bet<-reg$coefficients; nb=length(bet)
  Vb<- vcov(reg)/sigma(reg); #Vb=solve(t(x)%*%diag(W, nrow = length(W))%*%x)
  S = sum(W*(reg$residuals^2))
  set.seed(seed = seed); 
  if(sum(abs(W))<=nb){
    stop("Not enough degree of freedom for weighted error term variance")
  }
  sig<- 1/stats::rgamma(ndrs,shape = (sum(W)-nb)/2,scale = 2/S)
  obtbets<- function(sig){
    set.seed(seed)
    beta=MASS::mvrnorm(n=1,bet,sig*Vb)
    return(beta)
  }
  beta=sapply(sig,obtbets)
  val<- list(beta=t(beta),sig=sig)
  return(val)
}


#===========================================================================================#
#' Bayesian Normal Regression
#'
#' This function computes the bayesian normal regression with the option of using weights
#'
#' @param dat data frame with the first column being the dependent variable
#' @param useW a logical indicating if weights are to be used; default \code{useW=F}
#' @param xvec predictor vector with respect to which weights are taken
#' @param ndrs number of draws of coefficient and variance terms, default: \code{ndrs}=10000
#' @param seed the seed set for all random operations in the function, default \code{seed}=1
#' @return beta an ndrs x p matrix of coefficient draws, p is the number of OLS coefficients
#' @return sig a vector of variance (sigma squared) draws
#'
#'
#' @export

baysnreg<- function(dat,useW=F,xvec,ndrs=NULL,seed=NULL){
  if(useW){
    X = dat[,-1]
    #W = genmle::distnz(xvec,X)$kernW
  }else{
    W=NULL
  }
  if(is.null(ndrs)){
    ndrs=10000
  }
  if(is.null(seed)){
    seed=1
  }
  xx <- aclass(dat,ndrs=ndrs,seed=seed,W=W)
  rzlts<- simsigbet(xx)
  val<- list(beta=rzlts$beta,sigma=rzlts$sig)
}
#===========================================================================================#
#' Get predictive distribution
#'
#' This function obtains random draws from the predictive distribution
#'
#' @param baysobj object from \code{baysnreg} which has matrix of coefficient and \code{sigma} draws.
#' @param xvec the predictor vector for which the predictive distribution is sought
#' @return ydist a vector representing the predictive distribution.
#'
#' @export


predist<- function(baysobj,xvec){
  beta<- as.matrix(baysobj$beta);sigma=as.matrix(baysobj$sigma)
  yh<- beta%*%matrix(c(1,as.numeric(xvec)),ncol=1)
  ydist<- rnorm(length(sigma),mean=yh,sd=sqrt(sigma))
  return(ydist)
}
#===========================================================================================#
#' Bayesian Normal Regression
#'
#' This function computes the bayesian normal regression with the option of using weights
#'
#' @param dat data frame with the first column being the dependent variable
#' @param useW a logical indicating if weights are to be used; default \code{useW=F}
#' @param xvec predictor vector for which a predictive distribution is sought
#' @param ndrs number of draws of coefficient and variance terms, default: \code{ndrs}=10000
#' @param seed the seed set for all random operations in the function, default \code{seed}=1
#' @return beta an ndrs x p matrix of coefficient draws, p is the number of coefficients
#' @return ydist draws from the predictive distribution
#'
#' @examples
#' Ddat<- bayesprdopt::DAT_ALL[1:75,]; dat<- Ddat[-20,]; xvec=Ddat[20,-1];
#' dd<- predbreg(dat = dat,xvec = xvec); plot(density(dd$ydist)); plot(density(dd$beta[,3]))
#' ## A weighted version with non-default number of draws
#' dw<- predbreg(dat = dat,xvec = xvec,useW = TRUE,ndrs = 9000); plot(density(dw$ydist))
#'
#' @export

predbreg<- function(dat,xvec,useW=F,ndrs=NULL,seed=NULL){
baysobj<- baysnreg(dat,useW=useW,xvec=xvec,ndrs=ndrs,seed=seed)
  ydist<- predist(baysobj,xvec)
  val<- list(beta=baysobj$beta,sigma=baysobj$sigma,ydist=ydist)
  return(val)
}
#===========================================================================================#
#' @export

predbregy<- function(dat,xvec,useW=F,ndrs=NULL,seed=NULL){
  baysobj<- baysnreg(dat,useW=useW,xvec=xvec,ndrs=ndrs,seed=seed)
  ydist<- predist(baysobj,xvec)
  return(ydist)
}


#===========================================================================================#
#' @export
predbregj<- function(j,dat,xvec,grd,useW=F,ndrs=NULL,seed=NULL){
  idcls<- which(grd[j,]==1)
  y=dat[,1]; x<- dat[,-1][,idcls]; xvec=xvec[idcls]
  dat<- data.frame(y,x)
  return(predbregy(dat,xvec,useW=useW,ndrs=ndrs,seed=seed))
}
#===========================================================================================#
#' Bayesian predictive distribution with best subset selection
#'
#' This function computes predictive distributions using the linear regression model and
#' best subset selection
#'
#' @param dat data frame with the first column being the dependent variable
#' @param xvec predictor vector for which a predictive distribution is sought
#' @param mods a vector of indexes of models on the \code{grid} from \code{modslct} to use. The
#' default computes predictive distributions for all models.
#' @param sbset minimum number of covariates in each model.
#' Use negative to have maximum number of covariates in each model.
#' @param useW a logical indicating if weights are to be used; default \code{useW=F}
#' @param ndrs number of draws of coefficient and variance terms, default: \code{ndrs}=10000
#' @param seed the seed set for all random operations in the function, default \code{seed}=1
#' @return mat a matrix, with columns corresponding to predictive distributions
#'
#' @examples
#' Ddat<- bayesprdopt::DAT_ALL[1:75,]; dat<- Ddat[-20,]; xvec=Ddat[20,-1];
#' mdist<- predbregbs(dat,xvec,mods=NULL,sbset=9,useW=F,ndrs=NULL,seed=NULL)
#' plot(density(mdist[,4]))
#'
#' @export

predbregbs<- function(dat,xvec,mods=NULL,sbset=NULL,useW=F,ndrs=NULL,seed=NULL){
grd<- modslct(daten = dat[,-1],sbset = sbset)
if(is.null(mods)){
  mods<- 1:nrow(grd)
}
no_cores<-parallel::detectCores() - 1
if(length(mods)<=no_cores){
  mat<- sapply(mods,predbregj,dat=dat,xvec=xvec,grd=grd,useW=useW,ndrs=ndrs,seed=seed)
}else{
  c1<-parallel::makeCluster(no_cores, type = "PSOCK")
  mat<- parallel::parSapply(c1,mods,predbregj,dat=dat,xvec=xvec,grd=grd,useW=useW,ndrs=ndrs,seed=seed)
  parallel::stopCluster(c1)
}
return(mat)
}
#===========================================================================================#
#' Bayesian predictive distribution with supervised principal components
#'
#' @param dat data frame with the first column being the dependent variable
#' @param xvec predictor vector for which a predictive distribution is sought
#' @param mods a vector of indexes of models on the \code{grid} from \code{modslct} to use. The
#' default computes predictive distributions for all models.
#' @param sbset minimum number of covariates in each model.
#' Use negative to have maximum number of covariates in each model.
#' @param useW a logical indicating if weights are to be used; default \code{useW=F}
#' @param ndrs number of draws of coefficient and variance terms, default: \code{ndrs}=10000
#' @param seed the seed set for all random operations in the function, default \code{seed}=1
#' @return mat a matrix, with columns corresponding to predictive distributions
#'

#' @examples
#' Ddat<- bayesprdopt::DAT_ALL[1:75,]; dat<- Ddat[-20,]; xvec=Ddat[20,-1];
#' yspc<- predbregspc(dat,xvec); plot(density(yspc))
#' 
#' @export

predbregspc<- function(dat,xvec,useW=F,ndrs=NULL,seed=NULL,lev1=0.1,lev2=0.95,maxc=3, minft=6){
  xx<-supc(dat=dat,xvec = xvec,lev1=lev1,lev2=lev2,maxc=maxc, minft=minft);
  xvec=xx[1,]; x = xx[-1,]; dat=data.frame(dat[,1],x)
  ydist<- predbregy(dat,xvec,useW=useW,ndrs=ndrs,seed=seed)
  return(ydist)
}

#===========================================================================================#
#' @export
predspcj<- function(j,dat,useW = T){
  # this function can run for ARDL and Polynomial reg, both weighted and unweighted.
  xvec<- dat[j,-1]; ddat<- dat[1:(j-1),]
  yspcw<- predbregspc(ddat,xvec,useW = useW)
  return(yspcw)
}
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

#============================================
#' Get density values
#' 
#' Get density values for a continuous distribution whose random draws are the vector 
#' \code{vec}
#' 
#' @param vec continuous distribution whose random draws are elements of vec
#' @param val a vector of values in the support of vec for which to get pdfs
#' @param type the "probability" or "density" that is preferred as output
#' 
#' @return mval a vector of probability weights of each element of vec
#' 
#' @examples 
#' set.seed(40); v = rnorm(1000)**2; plot(density(v)); # chi-square distributed
#' # obtain pdfs at a vector of values x: 
#' x=c(0,2,4); dn<- getprobs(v,val=x); points(dn$x, dn$y, col = "red");dn
#' 
#' @export

getprobs<- function(vec,val=NULL,type="probability"){
  vec<- sort(vec)
  if(is.null(val)){
    val=vec
  }
  zz<- density(as.matrix(vec),n=1024)
  fsp =stats::splinefun(zz$x, zz$y)
  const<-integrate(fsp, min(zz$x), max(zz$x))
  if(type=="probability"){
    if(abs(const$value-1) < 1e-03){
      spdd<- stats::spline(zz$x,zz$y,xout = val)
    }else{
    fh1<- ks::kde(x=vec,binned = TRUE)
    d1<- ks::dkde(x=vec,fhat = fh1)
    spdd<- stats::spline(vec,d1,xout = val)
    }
    }else if(type=="density"){
    spdd<- stats::spline(zz$x,zz$y,xout = val)
    }else{
    stop(paste("Type ", type  ," not recognised. It is either \"probability\" or \"density\"",sep = ","))
    }
  mval = list(x=spdd$x,y=spdd$y)
  return(mval)
}
