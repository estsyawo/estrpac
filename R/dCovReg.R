#==========================================================================================>
#' MMD Regression
#'
#' \code{mmdreg.fit} runs a linear minimum mean dependence (MMD) regression.
#'
#' @param Y outcome variable
#' @param X matrix of covariates.
#' @param Z matrix of instruments. Defaults to \code{X}.
#' @param cl number of clusters to pass to \code{pbsapply()}. This is only advised in large samples.
#' @return an IV regression object which also contains coefficients, standard errors, etc. The
#' standard errors are computed based on a U-Statistics characterisation of the influence function
#' 
#' @importFrom stats dist
#' 
#' @examples 
#' ## Generate data and run MMD regression
#' n=200; set.seed(12); X = rnorm(n); er = rchisq(n,df=1)-1; Z=X; X=scale(abs(X))+er/sqrt(2)
#' Y=X+er
#' reg = mmdreg.fit(Y,X,Z) #run regression
#' ## MMD coefficients, standard errors, and t-statistics
#' reg$MMD_coefficients; reg$MMD_SE; reg$MMD_tstat
#' @export

mmdreg.fit = function(Y,X,Z=X,cl=NULL){
  YY = Y - mean(Y);XX=X=as.matrix(X); 
  for (k in 1:ncol(XX)){XX[,k]=X[,k]-mean(X[,k])}; n = length(Y)
  #in case Z is already a Euclidean distance matrix
  Z = as.matrix(Z)
  if(!(ncol(Z)==nrow(Z))){Mz = as.matrix(dist(Z))}else{Mz=Z}
  
  fn = function(i) apply(c(Mz[i,])*X,2,mean)
  if(is.null(cl)){Zhat = t(sapply(1:n,fn))}else{Zhat = t(pbapply::pbsapply(1:n,fn,cl=cl))}
  if(nrow(Zhat)!=n){Zhat=t(Zhat)} #in case row instead of column matrix
  Zhat = Zhat*n/(n-1) #U-statistic correction
  obj=ivreg::ivreg(YY~as.matrix(XX)|as.matrix(Zhat),x=TRUE)
  fn = function(i) sum(c(Mz[i,])*obj$residuals)/(n-1)
  if(is.null(cl)){Uhat = t(sapply(1:n,fn))}else{Uhat = t(pbapply::pbsapply(1:n,fn,cl=cl))}
  if(nrow(Uhat)!=n){Uhat=t(Uhat)} #in case row instead of column matrix
  A=crossprod(Zhat,XX)/n; B = crossprod(Zhat*obj$residuals + XX*c(Uhat))/n
  VC = solve(A,solve(A,B))/n
  obj$MMD_VC = VC; obj$MMD_SE = sqrt(diag(VC)); obj$MMD_coefficients=obj$coefficients[-1]
  obj$MMD_tstat=obj$MMD_coefficients/obj$MMD_SE; obj$MMD_Zm=Mz
  obj$MMD_Z=Z;class(obj)=c("MMD",class(obj))
  obj
}
#==========================================================================================>


#==========================================================================================>
#' Generic Linear Integrated Moment Regression
#'
#' \code{imlmreg.fit} runs a generic linear integrated moment regression allowing for different
#' kernels. 
#'
#' @param Y outcome variable
#' @param X matrix of covariates.
#' @param Z matrix of instruments
#' @param Kern type of kernel. See Details for available kernels
#' @param vctype type of sandwich covariance matrix (see \link[sandwich]{vcovHC})
#' 
#' @details The \eqn{(i,j)}'th elements of available kernel methods are
#' \describe{
#' \item{"Euclid"}{Euclidean distance between two vectors: ||Z_i-Z_j||} 
#' \item{"Gauss.W"}{The weighted Gaussian kernel: exp(-0.5(Z_i-Z_j)'V^{-1}(Z_i-Z_j)) where V is the variance of V} 
#' \item{"Gauss"}{The unweighted Gaussian kernel: exp(-||Z_i-Z_j||^2)}
#' \item{"DL"}{The kernel of Dominguez & Lobato 2004: \eqn{1/n\sum{l=1}^n I(Z_i\le Z_l)I(Z_j\le Z_l)}}
#' \item{"Esc6"}{The projected version of the DL in Escanciano 2006.}
#' \item{"WMD"}{The kernel used in Antoine & Lavergne 2014. See page 60 of paper.}
#' \item{"WMDF"}{The Fuller (1977)-like modification of the kernel in Antoine & Lavergne 2014. See page 64 of paper.}
#' }
#'  
#' @return an IV regression object which also contains coefficients, standard errors, etc.
#' 
#' @importFrom stats dist
#' @importFrom sandwich vcovHC
#' @importFrom ivreg ivreg
#' 
#' @examples 
#' ## Generate data and run MMD regression
#' n=200; set.seed(12); X = rnorm(n); er = (rchisq(n,df=1)-1)/sqrt(2); Z=X
#' X=scale(abs(X))+er/sqrt(2); Y=X+er
#' summary(imlmreg.fit(Y=Y,X=X,Z=Z))
#' summary(ivreg::ivreg(formula = Y ~ X | Z)) #compare to conventional IV regression
#' @export

imlmreg.fit = function(Y,X,Z,Kern="Euclid",vctype="HC3"){
  YY = Y - mean(Y);XX=as.matrix(X)
  for (k in 1:ncol(XX)){XX[,k]=XX[,k]-mean(XX[,k])}; n = length(Y)
  Z = as.matrix(Z) #in case Z is already a Euclidean distance matrix
  if(!(ncol(Z)==nrow(Z)&isSymmetric(Z))){
    Mz = Kern.fun(Z,Kern,XX,YY)
  }else{
      Mz=Z
      }
  Zhat = Mz%*%XX/(n-1)
  obj=ivreg::ivreg(YY~as.matrix(XX)|as.matrix(Zhat),x=TRUE); obj$Z=Z
  obj$vcovHC=sandwich::vcovHC(obj,type=vctype)
  obj$HC_Std.Err=sqrt(diag(obj$vcovHC))
  class(obj)=c(class(obj),"ICM",Kern)
  obj
}

#==========================================================================================>
#' ICM-IV Regression
#'
#' \code{imlmreg2.fit} runs a generic linear integrated moment regression allowing for different
#' kernels. This variant uses centred instruments in the meat of the sandwich matrix
#'
#' @param Y outcome variable
#' @param X matrix of covariates.
#' @param Z matrix of instruments
#' @param weights a vector of length \eqn{n} of weights for observations
#' @param Kern type of kernel. See Details for available kernels
#' @param vctype type of sandwich covariance matrix (see \link[sandwich]{vcovHC})
#' @param cluster vector of length \eqn{n} with cluster assignments of observations.
#' @param clus.est.type options are "A" and "B". "A" sets \eqn{K(Z_i,Z_j)=0} for \eqn{i,j} in the
#' same cluster while option "B" only does so for \eqn{i=j}.
#'
#' @details The \eqn{(i,j)}'th elements of available kernel methods are
#' \describe{
#' \item{"Euclid"}{Euclidean distance between two vectors: ||Z_i-Z_j||}
#' \item{"Gauss.W"}{The weighted Gaussian kernel: exp(-0.5(Z_i-Z_j)'V^{-1}(Z_i-Z_j)) where V is the variance of V}
#' \item{"Gauss"}{The unweighted Gaussian kernel: exp(-||Z_i-Z_j||^2)}
#' \item{"DL"}{The kernel of Dominguez & Lobato 2004: \eqn{1/n\sum{l=1}^n I(Z_i\le Z_l)I(Z_j\le Z_l)}}
#' \item{"Esc6"}{The projected version of the DL in Escanciano 2006.}
#' \item{"WMD"}{The kernel used in Antoine & Lavergne 2014. See page 60 of paper.}
#' \item{"WMDF"}{The Fuller (1977)-like modification of the kernel in Antoine & Lavergne 2014. See page 64 of paper.}
#' }
#'
#' @return an IV regression object which also contains coefficients, standard errors, etc.
#'
#' @importFrom stats dist
#' @importFrom sandwich vcovHC
#' @importFrom ivreg ivreg
#'
#' @examples
#' ## Generate data and run MMD regression
#' n=200; set.seed(12); X = rnorm(n); er = (rchisq(n,df=1)-1)/sqrt(2); Z=X
#' X=scale(abs(X))+er/sqrt(2); Y=X+er
#' summary(imlmreg2.fit(Y=Y,X=X,Z=Z))
#' summary(ivreg::ivreg(formula = Y ~ X | Z)) #compare to conventional IV regression
#' @export

imlmreg2.fit = function(Y,X,Z,weights=NULL,Kern="Euclid",vctype="HC0",
                        cluster=NULL,clus.est.type="A"){
  X=as.matrix(cbind(1,X))
  n = length(Y)
  Z = as.matrix(Z)
  
  if(!isSymmetric(Z)){
    Mz = Kern.fun(Z,Kern,X[,-1],Y)
  }else{
    Mz=Z
  }
  if(!is.null(cluster)){
    uclus<- unique(cluster); G<- length(uclus)
    ## tailor the Kernel matrix for cluster jackknifing ...
    if(clus.est.type=="A"){
      for(g in 1:G){
        idclus.g<- which(cluster==uclus[g])
        Mz[idclus.g,idclus.g]<- 0
      }
    }else if(clus.est.type=="B"){
      diag(Mz)<- 0.0
    }else{
      stop("Estimator type under clustered data not recognised.")
    }
    
  }#end if(!is.null(cluster))
  Zhat = Mz%*%X/(n-1)
  obj=ivreg::ivreg(Y~as.matrix(X[,-1])|as.matrix(Zhat[,-1]),x=TRUE,weights = weights)
  obj$Z=Z; obj$Mz = Mz
  if(is.null(cluster)){
    obj$vcovHC=sandwich::vcovHC(obj,type=vctype)
  }else{
    obj$vcovHC=sandwich::vcovCL(obj,cluster = cluster)
  }
  
  obj$HC_Std.Err=sqrt(diag(obj$vcovHC))
  class(obj)=c(class(obj),"ICM",Kern)
  obj
}
#==========================================================================================>


#==========================================================================================>
#' Generic k-class estimator
#'
#' \code{kClassIVreg.fit} runs a generic linear IV model of the k-Class
#'
#' @param Y outcome variable
#' @param X matrix of covariates.
#' @param Z matrix of instruments. Defaults to \code{X}.
#' @param method method of the k-class to implement. Defaults to "JIVE".
#' @param vctype type of sandwich covariance matrix (see \link[sandwich]{vcovHC})
#' @param cluster vector of length \eqn{n} with cluster assignments of observations.
#' @param weights a vector of length \eqn{n} of weights for observations
#' 
#' @return an IV regression object which also contains coefficients, standard errors, etc.
#' @details Available methods in the k-Class include
#' \describe{
#' \item{"JIVE"}{The Jackknife IV of Angrist et al. 1999} 
#' \item{"LIML"}{Limited Maximum Likelihood} 
#' \item{"HLIM"}{The Jackknife Limited Maximum Likelihood of Hausman et al. 2012}
#' \item{"HFUL"}{The heteroskedasticity robust version of the Fuller (1977) estimator}
#' }

#' @importFrom sandwich vcovHC
#' @importFrom ivreg ivreg
#' 
#' @examples 
#' ## Generate data and run MMD regression
#' n=200; set.seed(12); X = rnorm(n); er = rchisq(n,df=1)-1; Z=X 
#' Z=cbind(Z,Z^2,Z^3,Z^4);X=scale(abs(X))+er/sqrt(2); Y=X+er
#' summary(kClassIVreg.fit(Y=Y,X=X,Z=Z))
#' summary(ivreg::ivreg(formula = Y ~ X | Z)) #compare to conventional IV regression
#' @export

kClassIVreg.fit = function(Y,X,Z,method="JIVE",vctype="HC3",cluster=NULL,weights=NULL){
  n=length(Y);Z = as.matrix(Z) #in case Z is already an appropriate kernel matrix
  X=as.matrix(X) #the design matrix
  
  if(!(ncol(Z)==nrow(Z)&isSymmetric(Z))){
    Mz = Z%*%solve(crossprod(Z))%*%t(Z)
    if(method=="JIVE"){diag(Mz)=0 #Jackknife IV
    }else if(method=="LIML"){#Limited Information Maximum Likelihood
      Ystar = as.matrix(cbind(Y,1,X))
      diag(Mz)=diag(Mz)-min(Re(eigen(solve(crossprod(Ystar))%*%(crossprod(Ystar,Mz)%*%Ystar))$values))
    }else if(method=="HLIM"){
      Ystar = as.matrix(cbind(Y,1,X))
      diag(Mz)=-min(Re(eigen(solve(crossprod(Ystar))%*%(crossprod(Ystar,Mz)%*%Ystar))$values))
    }else if(method=="HFUL"){
      Ystar = as.matrix(cbind(Y,1,X))
      lam=min(Re(eigen(solve(crossprod(Ystar))%*%(crossprod(Ystar,Mz)%*%Ystar))$values))
      diag(Mz)=- (lam-(1-lam)/n)/(1-(1-lam)/n)
    }else{stop("The method ",method," is unavailable.\n")}
  }else{Mz=Z}
  
  Zhat = Mz%*%X/n
  obj=ivreg::ivreg(Y~X|as.matrix(Zhat),x=TRUE,y=TRUE,weights = weights); obj$Z=Z
  if(is.null(cluster)){
    obj$vcovHC=sandwich::vcovHC(obj,type=vctype)
  }else{
    obj$vcovHC=sandwich::vcovCL(obj,cluster = cluster)
  }
  obj$HC_Std.Err=sqrt(diag(obj$vcovHC))
  class(obj)=c(class(obj),"KClass",method)
  obj
}
#==========================================================================================>


#==========================================================================================>
#' Martingale Difference Divergence
#' 
#' \code{MDD} computes the Martingale Difference Divergence between univariate U and 
#' possibly multivariate Z. The MDD measures the mean dependence of U on Z.
#' 
#' @param U the univariate variable
#' @param Z the possibly multivariate conditioning variable
#' @return the MDD coefficient
#' 
#' @importFrom stats dist
#' 
#' @examples 
#' set.seed(12); X = rnorm(200); MDD(X,abs(X)); MDD(abs(X),X)
#' @export
MDD<- function(U,Z){
  n = length(U)
  U = U - mean(U); Z=as.matrix(Z)
  if(nrow(Z)!=n){Z=t(Z)}
  A = tcrossprod(U)
  if(is.matrix(Z)&isSymmetric(Z)&all(c(n==dim(Z)))){
    n=n #dummy operation
  }else{Z = as.matrix(dist(Z))}
  if(!all(n==c(dim(A),dim(Z)))){stop("mismatch in dimensions")}
  -sum(A*Z)/n^2
}
#==========================================================================================>
#===============================================================================#
#' Kernel matrix Esc6
#'
#' This a wrapper for a C function that computes the Kernel matrix of Escanciano 2006.
#'
#' @param Z n by pz matrix of instruments Z
#' @return n by n kernel matrix
#'
#' @examples
#' set.seed(12); Z = rnorm(5); Z=cbind(Z,abs(Z))
#' Kern.fun_Esc6(Z)
#' @useDynLib estrpac Kern_Esc
#' @export

Kern.fun_Esc6<- function(Z)
{
  Z = as.matrix(Z) # add intercept term
  n = as.integer(nrow(Z)); p = as.integer(ncol(Z))
  Z = as.double(Z); p = as.integer(p)
  Omg = as.double(matrix(0.0,n,n))
  ans=.C("Kern_Esc",Z,Omg=Omg,n,p)
  matrix(ans$Omg,n,n)
}
#===============================================================================#
#' Kernel matrix DL
#'
#' This is a wrapper for a C function that computes the Kernel matrix of Dominguez & Lobato 2004.
#'
#' @param Z n by pz matrix of instruments Z
#' @return n by n kernel matrix
#'
#' @examples
#' set.seed(12); Z = rnorm(5); Z=cbind(Z,abs(Z))
#' Kern.fun_DL(Z)
#' @useDynLib estrpac Kern_DL
#' @export

Kern.fun_DL<- function(Z)
{
  Z = as.matrix(Z) # add intercept term
  n = as.integer(nrow(Z)); p = as.integer(ncol(Z))
  Z = as.double(Z); p = as.integer(p)
  Omg = as.double(matrix(0.0,n,n))
  ans=.C("Kern_DL",Z,Omg=Omg,n,p)
  matrix(ans$Omg,n,n)
}
#===============================================================================#

#==========================================================================================>
#' Construction of n x n Kernel Matrices
#' 
#' \code{Kern.fun} computes the n x n matrix of kernels used in the integrated moment class
#' of estimators. 
#' 
#' @param Z n x p matrix of instrumental variables
#' @param Kern type of kernel desired. 
#' @param X matrix of endogenous covariates. Ought to be demeaned.
#' @param Y the outcome variable. Ought to be demeaned.
#' @details The \eqn{(i,j)}'th elements of available kernel methods are
#' \describe{
#' \item{"Euclid"}{Negative of the Euclidean distance between two vectors: -||Z_i-Z_j||} 
#' \item{"Gauss.W"}{The weighted Gaussian kernel: exp(-0.5(Z_i-Z_j)'V^{-1}(Z_i-Z_j)) where V is the variance of V} 
#' \item{"Gauss"}{The unweighted Gaussian kernel: exp(-||Z_i-Z_j||^2)}
#' \item{"DL"}{The kernel of Dominguez & Lobato 2004: \eqn{1/n\sum{l=1}^n I(Z_i\le Z_l)I(Z_j\le Z_l)}}
#' \item{"Esc6"}{The projected version of the DL in Escanciano 2006.}
#' \item{"WMD"}{The kernel used in Antoine & Lavergne 2014. See page 60 of paper.}
#' \item{"WMDF"}{WMD à la Fuller (1977) used in Antoine & Lavergne 2014. See page 64 of paper.}
#' }
#'  
#' @return the \eqn{n\times n} kernel matrix
#' @importFrom stats dist cov
#' 
#' @examples 
#' set.seed(12); X = rnorm(5); Z=cbind(X,abs(X)); Kern.fun(Z,"DL")
#' @export
Kern.fun = function(Z,Kern="Euclid",X=NULL,Y=NULL){
  Z=as.matrix(Z); n = nrow(Z)
  if(Kern=="Euclid"){Omg=-as.matrix(dist(Z))}
  else if(Kern=="Gauss.W"){
    Omg=exp(-0.5*as.matrix(dist(Z%*%expm::sqrtm(solve(cov(Z)))))^2)}
  else if(Kern=="Gauss"){Omg=exp(-0.5*as.matrix(dist(Z))^2);diag(Omg)=0.0}
  else if(Kern=="DL"){
    Omg=Kern.fun_DL(Z)
  }else if(Kern=="Esc6"){
    Omg=Kern.fun_Esc6(Z)
  }else if(Kern=="WMD"){
    if(is.null(X)&is.null(Y)){stop("This method requires that both X and Y be specified.")}
    Omg=(1/sqrt(2*pi))^ncol(Z)*exp(-0.5*as.matrix(dist(scale(Z)))^2)
    diag(Omg)=0;
    YX = as.matrix(cbind(Y,1,X))
    diag(Omg)=-min(Re(eigen(solve(crossprod(YX))%*%(crossprod(YX,Omg)%*%YX))$values))
  }else if(Kern=="WMDF"){
    if(is.null(X)&is.null(Y)){stop("This method requires that both X and Y be specified.")}
    Omg=(1/sqrt(2*pi))^ncol(Z)*exp(-0.5*as.matrix(dist(scale(Z)))^2)
    diag(Omg)=0;
    YX = as.matrix(cbind(Y,1,X)); lam=min(Re(eigen(solve(crossprod(YX))%*%(crossprod(YX,Omg)%*%YX))$values))
    diag(Omg)=-(lam-(1-lam)/n)/(1-(1-lam)/n)
  }
  else if(Kern=="Laplace"){Omg=exp(-as.matrix(dist(Z)))}
  #else if(Kern=="Cauchy"){Omg=exp(-0.5*as.matrix(dist(Z))^2)}
  else{stop("The method ",Kern," is not available.")}
  return(Omg)
}

#==========================================================================================>


#==========================================================================================>
#' Martingale Difference Divergence Matrix
#' 
#' \code{MDDM} computes the Martingale Difference Divergence Matrix between multivariate U and 
#' possibly multivariate Z of Lee & Shao 2018.
#' 
#' @param U the multivariate variable
#' @param Z the possibly multivariate conditioning variable. One can also directly supply
#' an n x n distance matrix
#' @return the MDDM matrix
#' 
#' @importFrom stats dist
#' 
#' @examples 
#' set.seed(12); U = rnorm(50); U=cbind(U,abs(U)); MDDM(U,2*pnorm(U[,1]))
#' MDDM(U,as.matrix(dist(2*pnorm(U[,1]))))
#' @export
MDDM<- function(U,Z){
  U = as.matrix(U); ncU=ncol(U); Z = as.matrix(Z)
  if(ncU<=1){stop("U must be multivariate.")}
  for (j in 1:ncU) {U[,j] = U[,j]-mean(U[,j])}
  if(ncol(Z)!=nrow(Z)){Z=as.matrix(dist(Z))}
  -crossprod(U,Z)%*%U/(nrow(U)^2)
}
#==========================================================================================>


#==========================================================================================>
#' Mammen's Wild Bootstrap Weights
#'
#' \code{wmat.mammen} generates Mammen's two-point weight for the wild bootstrap.
#'
#' @param n sample size
#' @param B number of wild bootstrap samples.
#' @param seed seed for reproducibility.
#' @param cluster vector of length n with cluster ids if cluster-robust wild-bootstrap is
#'  used
#' @return an n by B matrix of two-point weights
#' 
#' @export

wmat.mammen<- function(n,B=200,seed=NULL,cluster=NULL){ 
  if(!is.null(seed)) {set.seed(seed = seed)}
  p = (sqrt(5)+1)/(2*sqrt(5))
  a = -(sqrt(5)-1)/2; b = (sqrt(5)+1)/2
  if(!is.null(cluster)){uclus=unique(cluster);nC=length(uclus)
  wmat = matrix(NA,n,B)
  m=matrix(sample(c(a,b),nC*B,replace = TRUE,prob = c(p,1-p)),nrow = nC,ncol = B)
  for(j in 1:nC){idj = which(cluster==uclus[j])
  wmat[idj,]=m[rep(j,length(idj)),]
  }
  }else{
    wmat=matrix(sample(c(a,b),n*B,replace = TRUE,prob = c(p,1-p)),nrow = n,ncol = B)}
  wmat
}
#==========================================================================================>


#==========================================================================================>
#' MMD Specification Test
#'
#' \code{mmdlmspec.b_test} conducts the Su \& Zheng 2017 specification test on a linear model
#' estimated with the MMD estimator. It is based on the wild bootstrap.
#'
#' @param mmd.Obj is an MMD regression output from mmdreg.fit()
#' @param B number of wild bootstrap samples. Defaults to 199
#' @param wmat in case the n x B matrix of wild bootstrap weights are supplied by the user.
#' @param cl an integer to indicate number of child-processes in pbapply::pbsapply().
#' @param cluster vector of length n with cluster ids if cluster-robust wild-bootstrap is
#'  used
#' @return the p-value, test statistic, 
#' 
#' @examples 
#' ## Generate data and run MMD regression
#' n=100; set.seed(12); X = rnorm(n); er = rchisq(n,df=1)-1; Z=X; X=scale(abs(X))+er/sqrt(2)
#' Y=X+er; reg1 = mmdreg.fit(Y,X,Z); reg2 = mmdreg.fit(Y,X,X) #run regression
#' mmdlmspec.b_test(reg1); mmdlmspec.b_test(reg2) #test under the null and the alternative
#' ## MMD coefficients, standard errors, and t-statistics
#' 
#' @export

mmdlmspec.b_test<- function(mmd.Obj,B=199,wmat=NULL,cl=NULL,cluster=NULL){#returns p-value
  X = mmd.Obj$model[[2]]; k = ncol(X); n = nrow(X)
  Tn_SZo = n*MDD(mmd.Obj$residuals,Z=mmd.Obj$MMD_Zm)
  if(is.null(wmat)){wmat = wmat.mammen(n=n,B=B,seed = n,cluster=cluster)}
  
  fn<- function(j){
    Ystar = mmd.Obj$fitted.values + sqrt(n/(n-k-1))*mmd.Obj$residuals*wmat[,j]
    n*MDD(mmdreg.fit(Ystar,X=X,Z=mmd.Obj$MMD_Zm)$residuals,mmd.Obj$MMD_Zm)
  }
  if(is.null(cl)){TnSZ_boot=sapply(1:B,fn)}else{TnSZ_boot=pbapply::pbsapply(1:B,fn,cl=cl)}
  list(statistic=Tn_SZo,p.value=mean(TnSZ_boot>Tn_SZo))
}
#==========================================================================================>


#==========================================================================================>
#' MMD Relevance Test
#'
#' \code{mmdlmrelv.b_test} conducts the Su \& Zheng 2017 specification test on a linear model
#' estimated with the MMD estimator. It is based on the wild bootstrap.
#'
#' @param X1 an n x p1 matrix of endogenous covariates
#' @param X2 an n x p2 matrix of exogenous covariates
#' @param Z an n x pz matrix of instruments
#' @param B number of wild bootstrap samples. Defaults to 199
#' @param wmat in case the n x B matrix of wild bootstrap weights are supplied by the user.
#' @param cl an integer to indicate number of child-processes in pbapply::pbsapply().
#' @param cluster vector of length n with cluster ids if cluster-robust wild-bootstrap is
#'  used
#' @return the p-value, test statistic, 
#' 
#' @examples 
#' n=100; set.seed(12); X=rnorm(n); er=(rchisq(n,df=1)-1)/sqrt(2)
#' mmdlmrelv.b_test(X+er,abs(X),X^3) #multiple endogenous covariates
#' mmdlmrelv.b_test(cbind(X+er,X-(er)^3),Z=X^3) #multiple endogenous covariates
#' 
#' @export

mmdlmrelv.b_test<- function(X1,X2=NULL,Z,B=199,wmat=NULL,cl=NULL,cluster=NULL){
  X1=as.matrix(X1); p1 = ncol(X1)
  if(p1==1){
    if(!is.null(X2)){X2=as.matrix(X2) #exogenous covariates
    mmd.Obj = mmdreg.fit(X1,X2,Z)
    ans=mmdlmspec.b_test(mmd.Obj,B=B,wmat=wmat,cl=cl,cluster=cluster)
    }else{#no exogenous covariates
      ans=mddtest.boot(X1,Z,B=B,wmat=wmat,cl=cl,cluster=cluster)
    }#end if(!is.null(X2))
  }else{#multiple endogenous covariates
    X = as.matrix(cbind(X1,X2)) #full set of covariates
    Eobj=eigen(MDDM(X,Z))
    px=ncol(Eobj$vectors)
    #identify the covariate with max element in the eigen vector
    #idD=which.max(abs(c(Eobj$vectors[,px])))
    idD=which.max(abs(c(Eobj$vectors[1:p1,px])))
    mmd.Obj = mmdreg.fit(X[,idD],X[,-idD],Z)
    ans=mmdlmspec.b_test(mmd.Obj,B=B,wmat=wmat,cl=cl,cluster=cluster)
  }#end if(p1==1)
  return(ans)
}
#==========================================================================================>


#==========================================================================================>
#' ICM Specification Test for linear models
#'
#' \code{speclmb.test} conducts the Su \& Zheng 2017 specification test on a linear model
#' estimated with the MMD estimator. It is based on the wild bootstrap.
#'
#' @param reg.Obj is a regression output
#' @param Kern type of kernel desired
#' @param B number of wild bootstrap samples. Defaults to 199
#' @param wmat in case the n x B matrix of wild bootstrap weights are supplied by the user.
#' @param cl an integer to indicate number of child-processes in pbapply::pbsapply()
#' @param cluster vector of length n with cluster ids if cluster-robust wild-bootstrap is used
#'
#' @return the p-value, test statistic, 
#' 
#' @importFrom ivreg ivreg
#' @importFrom stats coef lm
#' 
#' @examples 
#' ## Generate data and run MMD regression
#' n=100; set.seed(12); X = rnorm(n); er = rnorm(n)
#' Y1=X+er/sqrt(1+X^2); Y2=X+abs(X)+er/sqrt(1+X^2)
#' speclmb.test(imlmreg2.fit(Y1,X,X)) #null
#' speclmb.test(imlmreg2.fit(Y2,X,X)) #alternative
#' 
#' @export

speclmb.test<- function(reg.Obj,Kern="Euclid",B=199,wmat=NULL,cl=NULL,cluster=NULL){
  if(length(class(reg.Obj))==1){ #for OLS and IV/2SLS
  if(class(reg.Obj)=="lm"){
    #check if the regression object contains y and x used.
    ifelse(is.null(reg.Obj$y)&is.null(reg.Obj$x),stop("reg.Obj needs to contain x and y."),1)
    regfn=function(Y,X,Z){lm(Y~X)}
    X=reg.Obj$x[,-1]; Y=reg.Obj$y
    Ker=Kern.fun(X,Kern = Kern)
  }else if(class(reg.Obj)=="ivreg"){
    #check if the regression object contains y, x, and z used.
    ifelse(is.null(reg.Obj$y)&is.null(reg.Obj$x),stop("reg.Obj needs to contain x, y, and z."),1)
    regfn=function(Y,X,Z){ivreg(Y~as.matrix(X)|as.matrix(Z))}
    X=reg.Obj$x$regressors[,-1]; Z=reg.Obj$x$instruments[,-1]
    Y=reg.Obj$y; Ker=Kern.fun(Z,Kern = Kern)
  }
  }else{ #for ICM and K-Class
    if(class(reg.Obj)[2]=="ICM"){
      #check if the regression object contains y, x, and z used.
      ifelse(is.null(reg.Obj$y)&is.null(reg.Obj$x$regressors) & is.null(reg.Obj$Z),stop("reg.Obj needs to contain x and y."),1)
      regfn=function(Y,X,Z){imlmreg2.fit(Y,X,Z=reg.Obj$Mz)}
      X=reg.Obj$x$regressors[,-1]; Z=reg.Obj$Z
      Y=reg.Obj$y; Ker=reg.Obj$Mz
    }else if(class(reg.Obj)[2]=="KClass"){
      #check if the regression object contains y, x, and z used.
      ifelse(is.null(reg.Obj$y)&is.null(reg.Obj$x$regressors)&is.null(reg.Obj$Z),stop("reg.Obj needs to contain x and y."),1)  
      regfn=function(Y,X,Z){kClassIVreg.fit(Y,X,Z,method = class(reg.Obj)[3])}
      X=reg.Obj$x$regressors[,-1]; Z=reg.Obj$Z
      Y=reg.Obj$y; Ker=Kern.fun(Z,Kern = Kern)
    }
  }#end if(length(.))
  
  k = length(coef(reg.Obj))-1; n = length(Y)
  Tn_o = c(t(reg.Obj$residuals)%*%Ker%*%reg.Obj$residuals/(n-1))
  
  if(is.null(wmat)){wmat = wmat.mammen(n=n,B=B,seed = n,cluster=cluster)}
  
  fn<- function(j){
    tryCatch(
      expr = {
        Ystar = reg.Obj$fitted.values + sqrt(n/(n-k-1))*reg.Obj$residuals*wmat[,j]
        reg.ObjStar=regfn(Ystar,X,Z)
        res=c(t(reg.ObjStar$residuals)%*%Ker%*%reg.ObjStar$residuals/(n-1))
        return(res)
      },
      error = function(e){return(NA)},
      warning = function(w){return(NULL)},
      finally = { }
    )#end tryCatch    
  }
  if(!is.na(Tn_o)){
    if(is.null(cl)){Tn_boot=sapply(1:B,fn)}else{Tn_boot=pbapply::pbsapply(1:B,fn,cl=cl)}
    ans=list(statistic=Tn_o,p.value=mean(Tn_boot[!is.na(Tn_boot)] > Tn_o))
  }else{
    ans=list(statistic=Tn_o,p.value=NA)
  }
  return(ans)
}
#==========================================================================================>


#==========================================================================================>
#' MDep objective function for non-linear models
#'
#' \code{mdep.nl} computes the objective function of the MDep estimator 
#' for a non-linear model
#'
#' @param theta parameter vector
#' @param U.fun user-written function of the disturbance function; evaluates to an n x 1 vector
#' @param Z n x pz matrix of instruments; should be supplied in case Z.m is null
#' @param Z.m D- or U-centred Euclidean matrix from instrument matrix Z - optional
#' @param sc adjusts the objective function - positive for minimisation, 
#' negative for maximisation
#' @return function value
#' 
#' @importFrom stats dist
#' 
#' @examples 
#' ## Generate data and run MMD regression
#' n=50; set.seed(12); X = rnorm(n); Y = X + (X^2-1)/sqrt(2)
#' U.fun=function(theta) Y^2-2*Y*(X*theta)-(X*theta)^2
#' mdep.nl(1,U.fun,X)
#' 
#' @export

mdep.nl=function(theta,U.fun,Z=NULL,Z.m=NULL,sc=1){
  if(is.null(Z.m)){Mat=energy::U_center(as.matrix(stats::dist(Z)));Z.m = c(Mat[lower.tri(Mat)])}
  sc*mean(c(dist(U.fun(theta)))*Z.m)
}
#==========================================================================================>


#==========================================================================================>
#' MMD objective function for non-linear models
#'
#' \code{mmd.nl} computes the objective function of the MMD estimator 
#' for a non-linear model
#'
#' @param theta parameter vector
#' @param U.fun user-written function of the disturbance function; evaluates to an n x 1 vector
#' @param Z n x pz matrix of instruments; should be supplied in case Z.m is null
#' @param Z.m Euclidean matrix from instrument matrix Z - optional
#' @param sc adjusts the objective function - positive for minimisation, 
#' negative for maximisation
#' @return function value
#' 
#' @importFrom stats dist
#' 
#' @examples 
#' ## Generate data and run MMD regression
#' n=50; set.seed(12); X = rnorm(n); Y = X + (X^2-1)/sqrt(2)
#' U.fun=function(theta) Y^2-2*Y*(X*theta)-(X*theta)^2
#' mmd.nl(1,U.fun,Z=X)
#' 
#' @export

mmd.nl=function(theta,U.fun,Z=NULL,Z.m=NULL,sc=1){
  if(is.null(Z.m)){Z.m = c(dist(Z))}
  U = U.fun(theta); U=U-mean(U); U.m = tcrossprod(U)
  -sc*mean(c(U.m[lower.tri(U.m)])*Z.m)
}
#==========================================================================================>


#==========================================================================================>
#' An MDD-based Test of Mean Independence by Wild Bootstrap
#'
#' \code{mddtest.boot} tests the mean independence of U conditional on V using the wild
#' bootstrap following Shao & Zhang 2014.
#'
#' @param U univariate variable
#' @param V possibly multivariate variable
#' @param B number of wild bootstrap samples. Defaults to 199
#' @param wmat in case the n x B matrix of wild bootstrap weights are supplied by the user.
#' @param cl an integer to indicate number of child-processes in pbapply::pbsapply().
#' @param cluster vector of length n with cluster ids if cluster-robust wild-bootstrap is
#'  used
#' @return test statistic and p-value
#' 
#' @examples 
#' set.seed(12); X = rnorm(200)
#' mddtest.boot(X,abs(X)); mddtest.boot(abs(X),X)
#' @export
#' 
mddtest.boot<- function(U,V,B=199,wmat=NULL,cl=NULL,cluster=NULL){
  U = U-mean(U); n = length(U); V = as.matrix(dist(V))
  mdd0 = n*MDD(U,V)
  if(is.null(wmat)){
    wmat=wmat.mammen(n,B=B,seed=0,cluster=cluster)
  }
  fn = function(j) {Ustar=U*wmat[,j]; n*MDD(Ustar,V)}
  if(is.null(cl)){Tns=sapply(1:B, fn)}else{Tns=pbapply::pbsapply(1:B, fn,cl=cl)}
  ans=list()
  ans$statistic = mdd0; ans$p.value=(1+sum(Tns>ans$statistic))/(1+B)
  ans
}
#==========================================================================================>

#==========================================================================================>
#' An ICM Test of Mean Independence by Wild Bootstrap
#'
#' \code{mindep.boot} tests the mean independence of U conditional on V using the wild
#' bootstrap
#'
#' @param U univariate variable
#' @param V possibly multivariate variable
#' @param B number of wild bootstrap samples. Defaults to 199
#' @param Kern type of kernel. See Details for available kernels
#' @param wmat in case the n x B matrix of wild bootstrap weights are supplied by the user.
#' @param cl an integer to indicate number of child-processes in pbapply::pbsapply().
#' @param cluster vector of length n with cluster ids if cluster-robust wild-bootstrap is
#'  used
#' @return test statistic and p-value
#' 
#' @examples 
#' set.seed(12); X = rnorm(200)
#' mindep.boot(X,abs(X)); mindep.boot(abs(X),X)
#' @export
#' 
mindep.boot<- function(U,V,B=199,Kern="Gauss",wmat=NULL,cl=NULL,cluster=NULL){
  U = U-mean(U); n = length(U)
  Ker = Kern.fun(V,Kern = Kern)
  Tn0 = c(n*t(U)%*%Ker%*%U/(n*(n-1)))
  if(is.null(wmat)){
    wmat=wmat.mammen(n,B=B,seed=0,cluster=cluster)
  }
  fn = function(j) {Ustar=U*wmat[,j]; n*t(Ustar)%*%Ker%*%Ustar/(n*(n-1))}
  if(is.null(cl)){Tns=sapply(1:B,fn)}else{Tns=pbapply::pbsapply(1:B,fn,cl=cl)}
  ans=list()
  ans$statistic = Tn0; ans$p.value=(1+sum(Tns>ans$statistic))/(1+B)
  ans
}
#==========================================================================================>


#==========================================================================================>
#' An MDD-based Test of Mean Independence by Permutation Bootstrap
#'
#' \code{mddtest.boot} tests the mean independence of U conditional on V using the wild
#' bootstrap following Shao & Zhang 2014.
#'
#' @param U univariate variable
#' @param V possibly multivariate variable
#' @param B number of wild bootstrap samples. Defaults to 199
#' @param PBmat in case the n x B matrix of permuted ids is supplied by the user.
#' @param cl an integer to indicate number of child-processes in pbapply::pbsapply().
#'  used
#' @return test statistic and p-value
#' 
#' @examples 
#' set.seed(12); X = rnorm(200)
#' mddtest.perm(X,abs(X)); mddtest.perm(abs(X),X)
#' @export
#' 
mddtest.perm<- function(U,V,B=199,PBmat=NULL,cl=NULL){
  U = U-mean(U); n = length(U); V = as.matrix(dist(V))
  mdd0 = n*MDD(U,V)
  if(is.null(PBmat)){
    PBmat=matrix(NA,n,B); for(j in 1:B){set.seed(j);PBmat[,j]=sample(1:n,n,replace = FALSE)}
  }
  fn = function(j) {n*MDD(U[PBmat[,j]],V)}
  if(is.null(cl)){Tns=sapply(1:B, fn)}else{Tns=pbapply::pbsapply(1:B, fn,cl=cl)}
  ans=list()
  ans$statistic = mdd0; ans$p.value=(1+sum(Tns>ans$statistic))/(1+B)
  ans
}
#==========================================================================================>


#==========================================================================================>
#' An MDD-based t-Test of Mean Independence
#'
#' \code{mdd.t_test} tests the mean independence of Y conditional on Z using a two-sided
#' t-test.
#'
#' @param Y univariate variable
#' @param Z possibly multivariate variable
#' @return t-statistic, standard error, and p-value of test
#' 
#' @importFrom stats pnorm
#' 
#' @examples 
#' set.seed(12); X = rnorm(200)
#' mdd.t_test(X,abs(X)); mdd.t_test(abs(X),X)
#' @export

mdd.t_test=function(Y,Z){
  Z = as.matrix(Z); Y=scale(Y); X = Y + apply(Z,1,mean); X=scale(X) #elicit an endogenous X
  MMDobj=mmdreg.fit(Y,X,Z)
  SE_MMD = sqrt(MMDobj$MMD_VC); obj=list()
  Tstat=MMDobj$coefficients[-1]/SE_MMD
  obj$statistic=c(Tstat); obj$stderr=c(SE_MMD)
  obj$p.value=2*(1-pnorm(abs(obj$statistic)))
  obj
}

#==========================================================================================>

#' #==========================================================================================>
#' #' An MDD-based t-Test of Mean Independence
#' #'
#' #' \code{mdd.t_test} tests the mean independence of U conditional on Z using a two-sided
#' #' t-test.
#' #'
#' #' @param U univariate variable
#' #' @param Z possibly multivariate variable
#' #' @return t-statistic, parameter estimate, standard error, and p-value
#' #' 
#' #' @importFrom stats pchisq
#' #' 
#' #' @examples 
#' #' set.seed(12); U = rnorm(200)
#' #' mdd.chi_test(U,abs(U)); mdd.chi_test(abs(U),U)
#' #' @export
#' 
#' mdd.chi_test=function(U,Z){
#'   Z = as.matrix(Z); U=scale(U); X = U + apply(Z,1,mean); X=scale(X) #elicit an endogenous X
#'   n = length(U); dn = n*(n-1)
#'   Om.Z=as.matrix(dist(Z))
#'   eta = c(crossprod(X,Om.Z)%*%X)
#'   Pi.n= (-Om.Z + (Om.Z%*%tcrossprod(X)%*%Om.Z)/eta)/(n*(n-1))
#'   diag(Pi.n)=0#center
#'   h.U = crossprod(U,Pi.n)
#'   dU = c(h.U)*c(U)
#'   del.n = sum(dU)
#'   ste.del.n=sqrt(sum((dU)^2))
#'   
#'   obj=list()
#'   
#'   obj$statistic=c(del.n/(ste.del.n))^2
#'   obj$delta.n=del.n; obj$stderr=c(ste.del.n)
#'   #obj$p.value=2*(1-pnorm(abs(obj$statistic)))
#'   obj$p.value=pchisq(obj$statistic,1,lower.tail = FALSE)
#'   obj
#' }
#' #==========================================================================================>
#' 
#' #==========================================================================================>
#' #' An MDD-based t-Test of Mean Independence
#' #'
#' #' \code{mdd.1t_test} tests the mean independence of U conditional on Z using a ... 
#' #' t-test.
#' #'
#' #' @param U univariate variable
#' #' @param Z possibly multivariate variable
#' #' @return t-statistic, parameter estimate, standard error, and p-value
#' #' 
#' #' @importFrom stats pchisq
#' #' 
#' #' @examples 
#' #' set.seed(12); U = rnorm(200)
#' #' mdd.1t_test(U,abs(U)); mdd.1t_test(abs(U),U)
#' #' @export
#' 
#' mdd.1t_test=function(U,Z){
#'   Z = as.matrix(Z); U=scale(U); X = U + apply(Z,1,mean); X=scale(X) #elicit an endogenous X
#'   n = length(U); dn = n*(n-1)
#'   Om.Z=as.matrix(dist(Z))
#'   eta = c(crossprod(X,Om.Z)%*%X)
#'   Dk = (Om.Z%*%tcrossprod(X)%*%Om.Z)/eta
#'   Pi.n= (-Om.Z + Dk)/dn
#'   diag(Dk)=0#center
#'   h.U = crossprod(U,Dk)
#'   dU = c(h.U)*c(U)/dn
#'   del.n = c(t(U)%*%Pi.n%*%U)
#'   ste.del.n=sqrt(sum(c(dU)^2))
#'   
#'   obj=list()
#'   obj$stat_1t=del.n/(ste.del.n)
#'   obj$stat_chi=obj$stat_1t^2
#'   
#'   obj$delta.n=del.n; obj$stderr=ste.del.n
#'   obj$p.value_1t=pnorm(-obj$stat_1t)
#'   obj$p.value_chi=pchisq(obj$stat_chi,1,lower.tail = FALSE)
#'   obj
#' }
#' #==========================================================================================>