#==========================================================================================>
#' MMD Regression
#'
#' \code{mmdreg} is runs a linear minimum mean dependence (MMD) regression.
#'
#' @param Y outcome variable
#' @param X matrix of covariates.
#' @param Z matrix of instruments. Defaults to \code{X}.
#' @param cl number of clusters to pass to \code{pbsapply()}. This is only advised in large samples.
#' @return an IV regression object which also contains coefficients, standard errors, etc.
#' 
#' @examples 
#' ## Generate data and run MMD regression
#' n=200; set.seed(12); X = rnorm(n); er = rchisq(n,df=1)-1; Z=X; X=scale(abs(X))+er/sqrt(2)
#' Y=X+er
#' reg = mmdreg(Y,X,Z) #run regression
#' ## MMD coefficients, standard errors, and t-statistics
#' reg$MMD_coefficients; reg$MMD_SE; reg$MMD_tstat
#' @export

mmdreg = function(Y,X,Z=X,cl=NULL){
  YY = Y - mean(Y);XX=X; for (k in 1:ncol(X)){XX[,k]=X[,k]-mean(X[,k])}; n = length(Y)
  #in case Z is already a Euclidean distance matrix
  Z = as.matrix(Z)
  if(!(ncol(Z)==nrow(Z))){Mz = as.matrix(dist(Z))}else{Mz=Z}
  
  fn = function(i) apply(c(Mz[i,])*X,2,mean)
  if(is.null(cl)){Zhat = t(sapply(1:n,fn))}else{Zhat = t(pbapply::pbsapply(1:n,fn,cl=cl))}
  if(nrow(Zhat)!=n){Zhat=t(Zhat)} #in case row instead of column matrix
  Zhat = Zhat*n/(n-1) #U-statistic correction
  obj=AER::ivreg(YY~as.matrix(XX)|as.matrix(Zhat))
  fn = function(i) sum(c(Mz[i,])*obj$residuals)/(n-1)
  if(is.null(cl)){Uhat = t(sapply(1:n,fn))}else{Uhat = t(pbapply::pbsapply(1:n,fn,cl=cl))}
  if(nrow(Uhat)!=n){Uhat=t(Uhat)} #in case row instead of column matrix
  A=crossprod(Zhat,XX)/n; B = crossprod(Zhat*obj$residuals + XX*c(Uhat))/n
  VC = solve(A,solve(A,B))/n
  obj$MMD_VC = VC; obj$MMD_SE = sqrt(diag(VC)); obj$MMD_coefficients=obj$coefficients[-1]
  obj$MMD_tstat=obj$MMD_coefficients/obj$MMD_SE; obj$MMD_Zm=Mz
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
#' @param Z the possibly conditioning multivariate variable
#' @return the MDD coefficient
#' 
#' @examples 
#' set.seed(12); X = rnorm(200); MDD(X,X^2/sqrt(2))
#' set.seed(12); Y = rchisq(200,1)/sqrt(2); MDD(Y,X)
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
  -sum(A*Z)/n
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
#' \code{mmd.lmspec_test} conducts the Su \& Zheng 2017 specification test a linear model
#' estimated with the MMD estimator.
#'
#' @param mmd.Obj is an MMD regression output from mmdreg()
#' @param B number of wild bootstrap samples. Defaults to 199
#' @param wmat in case the n x B matrix of wild bootstrap weights are supplied by the user.
#' @param cl an integer to indicate number of child-processes in pbapply::pbsapply().
#' @param cluster vector of length n with cluster ids if cluster-robust wild-bootstrap is
#'  used
#' @return the p-value, test statistic, 
#' 
#' @examples 
#' ## Generate data and run MMD regression
#' n=200; set.seed(12); X = rnorm(n); er = rchisq(n,df=1)-1; Z=X; X=scale(abs(X))+er/sqrt(2)
#' Y=X+er; reg1 = mmdreg(Y,X,Z); reg2 = mmdreg(Y,X,X) #run regression
#' mmd.lmspec_test(reg1); mmd.lmspec_test(reg2) #test under the null and the alternative
#' ## MMD coefficients, standard errors, and t-statistics
#' 
#' @export

mmd.lmspec_test<- function(mmd.Obj,B=199,wmat=NULL,cl=NULL,cluster=NULL){#returns p-value
  X = mmd.Obj$model[[2]]; k = ncol(X)
  Tn_SZo = MDD(mmd.Obj$residuals,Z=mmd.Obj$MMD_Zm)
  if(is.null(wmat)){wmat = bayesprdopt:::wmat.mammen(n=n,B=B,seed = n,cluster=cluster)}
  
  fn<- function(j){
    Ystar = mmd.Obj$fitted.values + sqrt(n/(n-k-1))*mmd.Obj$residuals*wmat[,j]
    MDD(mmdreg(Ystar,X=X,Z=mmd.Obj$MMD_Zm)$residuals,mmd.Obj$MMD_Zm)
  }
  if(is.null(cl)){TnSZ_boot=sapply(1:B,fn)}else{TnSZ_boot=pbapply::pbsapply(1:B,fn,cl=cl)}
  list(p_value=mean(TnSZ_boot>Tn_SZo),Tn=Tn_SZo)
}
#==========================================================================================>
#' MDep objective function for non-linear models
#'
#' \code{mdep.nl} computes the objective function of the MDep estimator 
#' for a non-linear model
#'
#' @param theta parameter vector
#' @param U.fun user-written function of the disturbance function; evaluates to an n x 1 vector
#' @param Z n x pz matrix of instruments; should be supplied in case Z.m is null
#' @param Z.m Euclidean matrix from instrument matrix Z - optional
#' @param sc -  adjusts the objective function - positive for minimisation, 
#' negative for maximisation
#' @return function value
#' 
#' @examples 
#' ## Generate data and run MMD regression
#' n=50; set.seed(12); X = rnorm(n); Y = X + (X^2-1)/sqrt(2)
#' U.fun=function(theta) Y^2-2*Y*(X*theta)-(X*theta)^2
#' mdep.nl(1,U.fun,X)
#' 
#' @export

mdep.nl=function(theta,U.fun,Z=NULL,Z.m=NULL,sc=1){
  if(is.null(Z.m)){Z.m = c(lower.tri(energy::U_center(as.matrix(dist(Z)))))}
  sc*mean(c(dist(U.fun(theta)))*Z.m)
}

# Function to compute the objective function of the MMD estimator for a non-linear model
# theta - parameter vector
# U.fun - a function of theta to evaluate the n x 1 vector of disturbances
# Z.m - Euclidean matrix - optional
# Z - n x pz matrix of instruments
# sc -  adjusts the objective function - positive for minimisation, negative for maximisation
# mmd.nl=function(theta,U.fun,Z.m=NULL,Z=NULL,sc=1){
#   if(is.null(Z.m)){Z.m = c(dist(Z)))}
#   U = U.fun(theta); U=U-mean(U); 
#   sc*mean(c(lower.tri(tcrossprod(U)))*Z.m)
# }

