#-------------------------------------------------------------------------------------$
# HARA utility function value for a single i'th draw of excess returns from pi1
#' @export
HARA = function(i,parm,rho,W0,R1,pi1){
  a=parm[1] #risky allocation
  th=parm[2] # cost parameter
  val = (rho/(1-rho))*((1-th)*W0*(R1 + a*pi1[i])/rho)^(1-rho)
  return(val) 
}
#-------------------------------------------------------------------------------------$
# Amount traded. 
# a - allocation (fractio) to risky asset, 
# W0- current dollar wealth, 
# A is current dollar holding in risky asset
#' @export
qq = function(a,W0,A) abs(a*W0 - A)
#-------------------------------------------------------------------------------------$

#-------------------------------------------------------------------------------------$
#' Expected HARA utility
#' 
#' \code{expHARA} expresses the expected form of the HARA utility for a two-asset portfolio, risky and riskless.
#' 
#' @param parm a parameter vector of allocation to risky asset and trading cost as a fraction of 
#' wealth
#' @param rho the coefficient of relative risk aversion
#' @param W0 current wealth
#' @param R1 gross return on the risk-less asset
#' @param pi1 a vector of random draws from the distribution of excess return on the risky asset
#' @return expected utility
#' 
#' @examples 
#' set.seed(30); pi1 = rnorm(5000,0.03954741,0.1868175);R1=1.0046; 
#' expHARA(parm = c(0.5,0),rho = 1.5,W0=500,R1=R1,pi1=pi1)
#' 
#' @export
expHARA= function(parm,rho,W0,R1,pi1) { 
  mean(sapply(1:length(pi1),HARA,parm = parm,rho = rho,W0=W0,R1=R1,pi1=pi1))
}
#-------------------------------------------------------------------------------------$
#' Certainty equivalent from HARA class
#' 
#' This function computes the certainty equivalent using \code{expHARA()}
#' @param parm a parameter vector of allocation to risky asset and trading cost as a fraction of 
#' wealth
#' @param rho the coefficient of relative risk aversion
#' @param W0 current wealth
#' @param R1 gross return on the risk-less asset
#' @param pi1 a vector of random draws from the distribution of excess return on the risky asset
#' @return expected utility
#' 
#' @examples 
#' set.seed(30); pi1 = rnorm(15000,0.03954741,0.1868175);R1=1.0046; 
#' ceHARA(parm = c(0.5,0),rho = 1.5,W0=500,R1=R1,pi1=pi1)
#' 
#' @export
ceHARA<- function(parm,rho,W0,R1,pi1){ 
  rho*((((1-rho)/rho)*expHARA(parm = parm,rho = rho,W0=W0,R1=R1,pi1=pi1))^(1/(1-rho)))
}
#-------------------------------------------------------------------------------------$

#-------------------------------------------------------------------------------------$
#' Distribution of Portfolio wealth
#' 
#' This function generates the distribution of wealth of a two-asset portfolio
#'
#' @param parm a parameter vector of allocation to risky asset and trading cost as a fraction of 
#' wealth 
#' @param W0 current wealth
#' @param R1 gross return on the risk-less asset
#' @param pi1 a vector of random draws from the distribution of excess return on the risky asset
#' @return wealth distribution
#' 
#' @examples 
#' set.seed(30); pi1 = rnorm(15000,0.03954741,0.1868175); R1=1.0046; 
#' w<- W1(parm = c(0.5,0),W0=500,R1=R1,pi1=pi1); plot(density(w))
#' 
#' @export
W1 = function(parm,W0,R1,pi1) {
  (1-parm[2])*W0*(R1 + parm[1]*pi1)
}
#-------------------------------------------------------------------------------------$


#-------------------------------------------------------------------------------------$
#' Trading cost function
#' 
#' Computes borrowing plus transaction costs using a quadratic cost function
#' 
#' @param parm a parameter vector of allocation to risky asset and trading cost as a fraction of 
#' wealth 
#' @param W0 current wealth
#' @param R1 gross return on the risk-less asset
#' @param A current dollar holding in the risky asset
#' @param cpars a two by three matrix whose columns are the quadratic cost parameters,
#' with the first column corresponding to the fixed cost. For proportional costs, set
#' the first and third columns to zero.
#' @return cost the dollar value of trading cost
#' 
#' @examples 
#' W0 = 500;A=250;R1=1.0046; cpars<- rbind(c(0.005, 0.0025, 0.00025), c(0.001, 0.0005, 0.00005))
#' Costfn(parm = c(0.5,0),W0=W0,R1=R1,A=A,cpars = cpars)
#'
#' @export

Costfn = function(parm,W0,R1,A,cpars){
  alf<- cpars[,1]; bet<- cpars[,2]; gam<- cpars[,3]
  a<-parm[1]; amt<- qq(a,W0,A)
  cost = sum(alf) + amt*(.5*(abs(a)+abs(1-a)-1)*R1 + sum(bet)+amt*sum(gam) )
  return(cost)
}

#Cost constraint
#' @export
CostConst = function(parm,W0,R1,A,cpars){
  th = parm[2]
  val= Costfn(parm,R1,W0,A,cpars) - th*W0 
  return(-val)
}

#-------------------------------------------------------------------------------------$
#' Generate distribution of decision parameters
#' 
#' This function does constrained HARA utility maximisation and outputs optimal parameter values and
#' the joint distribution of parameters.
#'
#' @param start suggested starting values
#' @param rho coefficient of relative risk aversion
#' @param W0 current wealth
#' @param R1 gross return on the risk-less asset
#' @param pi1 a vector of random draws from the distribution of excess return on the risky asset
#' @param A current dollar holding in the risky asset
#' @param cpars a two by three matrix whose columns are the quadratic cost parameters,
#' with the first column corresponding to the fixed cost. For proportional costs, set
#' the first and third columns to zero.
#' @param ub upper bound on borrowing ; a fraction of wealth
#' @param uth upper bound on cost fraction parm[2]
#' @param scale a value multiplied by \code{propob$var} in order to adjust the proposal distribution.
#' The default is \code{1.1} but the user can adjust it until a satisfactory acceptance rate is obtained.
#' @param iter number of random draws desired (default: 15000)
#' @param burn burn-in period for the Random Walk MH algorithm (default: 1000)
#' @return par parameter values from constrained optimisation
#' @return var variance-covariance matrix of parameters
#' @return pardraws matrix of parameter draws
#' @return postvals distribution of utility obtained from \code{pardraws}
#' @return AcceptRatio the acceptance ratio
#' 
#' @examples 
#' ## dj<- optu(start = c(0.5,0),rho = rho,W0=W0,R1=R1,pi1 = pi1,A=A,cpars = cpars,ub=ub,uth = uth)
#' 
#' @export

optu<- function(start,rho,W0,R1,pi1,A,cpars,ub,uth,scale = 1.1, iter = 15000, burn = 1000){
  #Objective function
  expHARAfn<- function(parm) expHARA(parm,rho = rho,W0=W0,R1=R1,pi1=pi1)
  
  #Equality Cost constraint
  hheq=function(parm) { CostConst(parm,R1,W0,A,cpars)}
  
  #Inequality constraints on parameters
  hhin = function(parm) { #each has to be non-negative
    h = rep(NA,4)
    h[1] = uth - parm[2]
    h[2] = parm[2]
    h[3] = parm[1] + ub
    h[4] = - parm[1] + (1+ub)
    h
  }
  #Constrained optimisation
  uumax<- auglag(par=start, fn = expHARAfn, hin = hhin, heq = hheq,
                 control.outer = list(trace=FALSE,method="nlminb",fnscale=-1))
  #posterior (with substituted multipliers)
  utpost<- function(parm) {
    exp(expHARAfn(parm) + sum(uumax$lambda * c(hhin(parm),hheq(parm)))) 
  }
  H = solve(-uumax$hessian)
  gr = numDeriv::grad(utpost,uumax$par); gr = as.matrix(gr)
  #laplace approximation of objective function
  propob<- list(mode=uumax$par,var=H%*%(gr%*%t(gr))%*%H)
  if(any(eigen(propob$var)$values<0)){
    propob$var<- Matrix::nearPD(propob$var)$mat #force var-cov to positive definite if it's not
  }
  
  #Obtain parameter draws and function value distribution
  MHdraw<- IndepMHgen(start = NULL, posterior = utpost, propob = propob,
                     const = hhin, scale = scale, iter = iter, burn = burn)
  #values to return
  val<-list(par=uumax$par,var=var,pardraws=MHdraw$Matpram,postvals=MHdraw$postvals,
            AcceptRatio=MHdraw$AcceptRatio)
  return(val)
}









#-------------------------------------------------------------------------------------$
# Utility distribution
#-------------------------------------------------------------------------------------$


#-------------------------------------------------------------------------------------$
# HARA utility function value for a single i'th draw of excess returns from pi1, with equality
# cost fraction constraint substituted into the objective function
# Second element of parm is a dummy input. The actual value i.e. cost fraction is determined
# inside the fuction
# example: HAReq(i=2,parm=0.5,rho=1.5,W0=W0,R1=R1,pi1=pi1,A=A,cpars=cpars)
#' @export
HAReq <- function(i,parm,rho,W0,R1,pi1,A,cpars){
  th = Costfn(parm,W0,R1,A,cpars)/W0 #cost as function of parm
  a=parm[1] #risky allocation
  #th=parm[2] # cost parameter
  val = (rho/(1-rho))*((1-th)*W0*(R1 + a*pi1[i])/rho)^(1-rho)
  return(val) 
}

# vector of utility values from the distribution of utility
# Example: Uv = disHAReq(parm=0.5,rho=1.5,W0=W0,R1=R1,pi1=pi1,A=A,cpars=cpars)
disHAReq= function(parm,rho,W0,R1,pi1,A,cpars) { 
  Np = length(pi1)
  Uv=sapply(1:Np,HAReq,parm = parm,rho = rho,W0=W0,R1=R1,pi1=pi1,A=A,cpars=cpars)
  return(Uv)
}

# log-likelihood of utility
# example: llkHAReq(parm=0.5,rho=1.5,W0=W0,R1=R1,pi1=pi1,A=A,cpars=cpars)
llkHAReq<- function(parm,rho,W0,R1,pi1,A,cpars){
  dist<- disHAReq(parm=parm,rho=rho,W0=W0,R1=R1,pi1=pi1,A=A,cpars=cpars)
  dist = dist[!is.na(dist)]
  prbu <- getprobs(dist); prbu=prbu[which(prbu>0)]; dist=dist[which(prbu>0)]
  Np = length(prbu)
  val = sum(dist) + sum(log(prbu)) - Np*log(sum(exp(dist*prbu)))
  return(val)
}
#lku = function(parm) llkHAReq(parm,rho=1.5,W0=W0,R1=R1,pi1=pi1,A=A,cpars=cpars)
# lku(parm=0.2)

# example: parm=0.5;W0=W0;R1=R1;A=A;cpars=cpars;uth=5/100;ub=1.5; hhin(parm)
hhin = function(parm) { #each has to be non-negative
  h = rep(NA,3)
  h[1] = uth - Costfn(parm,W0,R1,A,cpars)/W0
  #h[2] = parm[2]
  h[2] = parm[1] + ub
  h[3] = - parm[1] + (1+ub)
  h
}

# obtain mode of the posterior distribution
# uumax<-alabama::auglag(par=parm, fn = lku, hin = hhin,
#                control.outer = list(trace=FALSE,method="nlminb",fnscale=-1))

#cfrac=  Costfn(parm = uumax$par,W0=W0,R1=R1,A=A,cpars = cpars)/W0 # cost fraction
# normal approx:eg: zg=rnorm(30,uumax$par,sd = sqrt(1/uumax$hessian))


# propob<- list(mode=uumax$par,var=1/uumax$hessian)
# 
# #Obtain parameter draws and function value distribution
# scale = 2; burn = 1000; iter = 11000
# MHdraw<- IndepMHgen(start = NULL, posterior = lku, propob = propob,
#                     const = hhin, scale = 2, iter = iter, burn = burn)

