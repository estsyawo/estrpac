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
  sig<-invgamma::rinvgamma(ndrs,shape = (nr-nb)/2,scale = S/2)

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
  y = dat[,1]; x = as.matrix(dat[,-1]); nr=nrow(x)
  reg<-lm(y~as.matrix(x),weights = W)
  bet<-reg$coefficients; nb=length(bet)
  Vb<- vcov(reg)/sigma(reg);
  S = sum(reg$residuals^2)
  set.seed(seed = seed)
  sig<-invgamma::rinvgamma(ndrs,shape = (nr-nb)/2,scale = S/2)
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
#' @importFrom genmle distnz
#'
#' @export

baysnreg<- function(dat,useW=F,xvec,ndrs=NULL,seed=NULL){
  if(useW){
    X = dat[,-1]
    W = genmle::distnz(xvec,X)$kernW
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
#' Ddat<- prdtrade::DAT_ALL[1:75,]; dat<- Ddat[-20,]; xvec=Ddat[20,-1];
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
#' Ddat<- prdtrade::DAT_ALL[1:75,]; dat<- Ddat[-20,]; xvec=Ddat[20,-1];
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
#' dat<- prdtrade::DAT_ALL[1:75,]; dat<- Ddat[-20,]; xvec=Ddat[20,-1];
#' yspc<- predbregspc(dat,xvec); plot(density(yspc))

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
#' \code{IndepMH} computes random draws of parameters using a normal proposal distribution.
#' This function implements a generic form of \code{IndepMH} from the package \code{bayesdistreg}
#'
#' @param start starting values of parameters for the MH algorithm.
#' It is automatically generated from the normal proposal distribution but the user can also specify.
#' @param posterior the posterior distribution function. 
#' @param ... additional arguments to the posterior function
#' Should take parameter input of the same length as \code{start} or \code{propob$mode}
#' @param propob a list of mode and variance-covariance matrix of the normal proposal distribution. 
#' Save list as propob=list(mode=mode,var=variance-covariance)
#' @param scale a value multiplied by \code{propob$var} in order to adjust the proposal distribution.
#' The default is \code{1.5} but the user can adjust it until a satisfactory acceptance rate is obtained.
#' @param iter number of random draws desired (default: 15000)
#' @param burn burn-in period for the Random Walk MH algorithm (default: 1000)
#' @return val a list of matrix of draws pardraws and the acceptance rate
#'
#' @export

IndepMH<- function(start=NULL,posterior=NULL,...,propob=NULL,scale=1.5,iter=15000,burn=1000){
  varprop = scale*propob$var
  npar = length(propob$mode)
  Mat = array(0, c(iter, npar))
  if(is.null(start)){
    start = MASS::mvrnorm(n=1,propob$mode,varprop)
  }
  Mat[1,] = start; AccptRate<-0
  for(i in 2:iter){
    start= Mat[i-1,]
    prop = MASS::mvrnorm(n=1,propob$mode,varprop)#make a draw from proposal dist
    lpa = posterior(prop,...); lpb = posterior(start,...)
    accprob = exp(lpa-lpb)
    # the other part cancels out because the normal distribution is symmetric
    if(stats::runif(1)< accprob){
      Mat[i,]=prop
      AccptRate<- AccptRate +1
    }else{
      Mat[i,]=start
    }
  }
  cat("IndepMH algorithm successful\n")
  val = list(Matpram=Mat[-c(1:burn),],AcceptanceRate = AccptRate/iter)
  return(val)
}

