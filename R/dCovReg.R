#===========================================================================================#
#' MMD Regression
#'
#' \code{mddreg} is runs a linear minimum mean dependence (MMD) regression.
#'
#' @param Y outcome variable
#' @param X matrix of covariates.
#' @param Z matrix of instruments. Defaults to \code{X}.
#' @param cl number of clusters to pass to \code{pbsapply()}. This is only advised in large samples.
#' @return an IV regression object which also contains coefficients, standard errors, etc.
#' 
#' @examples 
#' n=200; set.seed(n); X = rnorm(n); er = rchisq(n,df=1)-1; Z=X; X=scale(abs(X))+er/sqrt(2)
#' Y=X+er
#' reg = mddreg(Y,X,Z) #run regression
#' ## MMD coefficients, standard errors, and t-statistics
#' reg$MMD_coefficients; reg$MMD_SE; reg$MMD_tstat
#' @export

mddreg = function(Y,X,Z=X,cl=NULL){
  YY = Y - mean(Y);XX=X; for (k in 1:ncol(X)){XX[,k]=X[,k]-mean(X[,k])}; n = length(Y)
  Mz = as.matrix(dist(Z))
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
  obj$MMD_tstat=obj$MMD_coefficients/obj$MMD_SE
  obj
}
