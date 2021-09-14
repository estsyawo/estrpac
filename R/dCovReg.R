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
  obj=AER::ivreg(YY~as.matrix(XX)|as.matrix(Zhat),x=TRUE)
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
#' @param D matrix of endogenous covariates
#' @param X matrix of exogenous covariates.
#' @param Z matrix of instruments
#' @param Kern type of kernel. See \code{\link{Kern.fun}} for available kernels
#' @return an IV regression object which also contains coefficients, standard errors, etc.
#' 
#' @importFrom stats dist
#' @importFrom sandwich vcovHC
#' @importFrom AER ivreg
#' 
#' @examples 
#' ## Generate data and run MMD regression
#' n=200; set.seed(12); D = rnorm(n); er = rchisq(n,df=1)-1; Z=D
#' D=scale(abs(D))+er/sqrt(2); Y=D+er
#' summary(imlmreg.fit(Y=Y,D=D,Z=Z))
#' summary(AER::ivreg(formula = Y ~ D | Z)) #compare to conventional IV regression
#' @export

imlmreg.fit = function(Y,D,X=NULL,Z,Kern="Euclid"){
  YY = Y - mean(Y);XX=as.matrix(cbind(D,X))
  for (k in 1:ncol(XX)){XX[,k]=XX[,k]-mean(XX[,k])}; n = length(Y)
  Z = as.matrix(Z) #in case Z is already a Euclidean distance matrix
  if(!(ncol(Z)==nrow(Z)&isSymmetric(Z))){Mz = Kern.fun(Z,Kern,D-mean(D),YY)}else{Mz=Z}
  Zhat = Mz%*%XX/n
  obj=AER::ivreg(YY~as.matrix(XX)|as.matrix(Zhat))
  obj$vcovHC=sandwich::vcovHC(obj)
  obj$HC_Std.Err=sqrt(diag(obj$vcovHC))
  obj
}
#==========================================================================================>
#' Generic k-class estimator
#'
#' \code{kClassIVreg.fit} runs a generic linear IV model of the k-Class
#'
#' @param Y outcome variable
#' @param D matrix of endogenous covariates
#' @param X matrix of exogenous covariates.
#' @param Z matrix of instruments. Defaults to \code{X}.
#' @param method method of the k-class to implement. Defaults to "JIVE".
#' @return an IV regression object which also contains coefficients, standard errors, etc.
#' @details Available methods in the k-Class include
#' \describe{
#' \item{"JIVE"}{The Jackknife IV of Angrist et al. 1999} 
#' \item{"LIML"}{Limited Maximum Likelihood} 
#' \item{"HLIM"}{The Jackknife Limited Maximum Likelihood of Hausman et al. 2012}
#' }

#' @importFrom sandwich vcovHC
#' @importFrom AER ivreg
#' @importFrom expm sqrtm
#' 
#' @examples 
#' ## Generate data and run MMD regression
#' n=200; set.seed(12); D = rnorm(n); er = rchisq(n,df=1)-1; Z=D 
#' Z=cbind(Z,Z^2,Z^3,Z^4);D=scale(abs(D))+er/sqrt(2); Y=D+er
#' summary(kClassIVreg.fit(Y=Y,D=D,Z=Z))
#' summary(AER::ivreg(formula = Y ~ D | Z)) #compare to conventional IV regression
#' @export

kClassIVreg.fit = function(Y,D,X=NULL,Z,method="JIVE"){
  n=length(Y);Z = as.matrix(Z) #in case Z is already an appropriate kernel matrix
  X=as.matrix(cbind(D,X)) #the design matrix
  
  if(!(ncol(Z)==nrow(Z)&isSymmetric(Z))){
    Mz = Z%*%solve(crossprod(Z))%*%t(Z)
    if(method=="JIVE"){diag(Mz)=0 #Jackknife IV
    }else if(method=="LIML"){#Limited Information Maximum Likelihood
      Ystar = as.matrix(cbind(Y,1,D)); YSinv=sqrtm(solve(crossprod(Ystar)))
      diag(Mz)=diag(Mz)-min(eigen(YSinv%*%(crossprod(Ystar,Mz)%*%Ystar)%*%YSinv)$values)
    }else if(method=="HLIM"){
      Ystar = as.matrix(cbind(Y,1,D)); YSinv=sqrtm(solve(crossprod(Ystar)))
      diag(Mz)=-min(eigen(YSinv%*%(crossprod(Ystar,Mz)%*%Ystar)%*%YSinv)$values)
    }else{stop("The method ",method," is unavailable.\n")}
  }else{Mz=Z}
  
  Zhat = Mz%*%X/n
  obj=AER::ivreg(Y~X|as.matrix(Zhat))
  obj$vcovHC=sandwich::vcovHC(obj)
  obj$HC_Std.Err=sqrt(diag(obj$vcovHC))
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

#==========================================================================================>
#' Construction of n x n Kernel Matrices
#' 
#' \code{Kern.fun} computes the n x n matrix of kernels used in the integrated moment class
#' of estimators. 
#' 
#' @param Z n x p matrix of instrumental variables
#' @param Kern type of kernel desired. 
#' @param D matrix of endogenous covariates. Ought to be demeaned.
#' @param Y the outcome variable. Ought to be demeaned.
#' @details The \eqn{(i,j)}'th elements of available kernel methods are
#' \describe{
#' \item{"Euclid"}{Euclidean distance between two vectors: ||Z_i-Z_j||} 
#' \item{"Gauss.W"}{The weighted Gaussian kernel: exp(-0.5(Z_i-Z_j)'V^{-1}(Z_i-Z_j)) where V is the variance of V} 
#' \item{"Gauss"}{The unweighted Gaussian kernel: exp(-||Z_i-Z_j||^2)}
#' \item{"DL"}{The kernel of Dominguez & Lobato 2004: \eqn{1/n\sum{l=1}^n I(Z_i\le Z_l)I(Z_j\le Z_l)}}
#' \item{"WMD"}{The kernel used in Antoine & Lavergne 2014. See page 60 of paper.}
#' }
#'  
#' @return the \eqn{n\times n} kernel matrix
#' @importFrom stats dist cov
#' 
#' @examples 
#' set.seed(12); X = rnorm(5); Z=cbind(X,abs(X)); Kern.fun(Z,"DL")
#' @export
Kern.fun = function(Z,Kern="Euclid",D=NULL,Y=NULL){
  Z=as.matrix(Z); n = nrow(Z)
  if(Kern=="Euclid"){Omg=as.matrix(dist(Z))}
  else if(Kern=="Gauss.W"){
    Omg=exp(-0.5*as.matrix(dist(Z%*%expm::sqrtm(solve(cov(Z)))))^2)}
  else if(Kern=="Gauss"){Omg=exp(-0.5*as.matrix(dist(Z))^2)}
  else if(Kern=="DL"){
    Omg=Om=matrix(NA,n,n)
    fn=function(i){ 
      sub.fn=function(j){all(Z[i,]<=Z[j,])*1}; sub.fn=Vectorize(sub.fn)
      sub.fn(1:nrow(Z))}#end function fn
    for(i in 1:n){Om[i,]=fn(i)}
    for(i in 1:n){for(j in 1:i){Omg[i,j]=mean(Om[i,]*Om[,j]);Omg[j,i]=Omg[i,j]}}
  }else if(Kern=="WMD"){
    if(is.null(D)&is.null(Y)){stop("This method requires that both D and Y be specified.")}
    Omg=(1/sqrt(2*pi))^ncol(Z)*exp(-0.5*as.matrix(dist(scale(Z)))^2)
    diag(Omg)=0;
    YX = as.matrix(cbind(Y,D)); YSinv=sqrtm(solve(crossprod(YX)))
    diag(Omg)=-min(eigen(YSinv%*%(crossprod(YX,Omg)%*%YX)%*%YSinv)$values)
  }
  # else if(Kern=="HMMD"){
  #   if(is.null(D)&is.null(Y)){stop("This method requires that both X and Y be specified.")}
  #   Omg=as.matrix(dist(Z))
  #   YX = as.matrix(cbind(Y,D)); YSinv=sqrtm(solve(crossprod(YX)))
  #   diag(Omg)=-min(eigen(YSinv%*%(crossprod(YX,Omg)%*%YX)%*%YSinv)$values)
  # }
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
#' n=100; set.seed(12); X = rnorm(n); er = (rchisq(n,df=1)-1)/sqrt(2)
#' X1=(X+er)/sqrt(2);mmdlmrelv.b_test(X1,abs(X),X^3);mmdlmrelv.b_test(X1,X,X^3)
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
    #...
  }#end if(p1==1)
  
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