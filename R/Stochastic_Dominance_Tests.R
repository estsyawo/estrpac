#----------------------------------------------------
#' @export
classifSD<- function(X,Y,order,type,B){# a generic class creator
    val = structure(list(X=X,Y=Y,B=B), class=paste(c(type,toString(order)),collapse = ""))
    return(val)
}
#----------------------------------------------------
#' @export
SD_semipar<- function(X,Y){
  lX = length(X); lY=length(Y)
  N = min(lX,lY)#shorter vector length
  if(N>=10000){
  k = 100
  }else if(N<36){
    stop("Not enough observations in shorter vector")
  }else{
    k=floor(sqrt(N))
  }
  Dv = quantile(c(X,Y),probs = seq(0,1,length.out = (k+1))) #k intervals over support of X,Y
  D = abs(Dv[2:(k+1)]-Dv[1:k]);  D=unname(D)
  nX = rep(0,k); nY=rep(0,k)
  If = matrix(0,nrow = k,ncol = k); If[lower.tri(If,diag = T)]<- 1
  
  IF =oMat= matrix(0,nrow = k,ncol = k)
  IF[k,k] = D[k]
  for(j in 1:(k-1)){
    IF[which(c(1:k)>j),j]=D[j]+D[j+1]
    IF[j,j]=D[j]
    nX[j] = length(which(X>Dv[j] & X<=Dv[j+1]))
    nY[j] = length(which(Y>Dv[j] & Y<=Dv[j+1]))
  }
  nX[k]=length(which(X>Dv[k]));nY[k]=length(which(Y>Dv[k]))
  pX = nX/lX; pY = nY/lY; p = (nX + nY)/(lX+lY)
  for(i in 1:k){
    for(j in 1:k){
      if(i==j){
        oMat[i,j]= p[i]*(1-p[j])  
      }else{
        oMat[i,j]=-p[i]*p[j]
      }
    }
  }
  IF=0.5*IF
  moMat = ((lX+lY)/(lX*lY))*oMat
  v = pX - pY; v = matrix(v,ncol = 1)
  #First-Order Stochastic dominance:
  vf = If%*%v; vf = matrix(vf,ncol = 1)
  M_inv1 = MASS::ginv(t(If)%*%moMat%*%If)
  Tf = t(vf)%*%M_inv1%*%vf; 
  pval1=pchisq(Tf,df=(k-1)); 
  Alternative1=ifelse(any(vf>1e-06) & any(vf< -1e-06),"Indeterminate","First Order Stochastic Dominance") 
  First_Order = list(Statistic = Tf,P_value = pval1,Alternative=Alternative1)
  
  #Second-Order Stochastic dominance:
  vF = IF%*%If%*%v; vF=matrix(vF,ncol = 1)
  M_inv2 = MASS::ginv(t(IF%*%If)%*%moMat%*%(IF%*%If))
  TF = t(vF)%*%M_inv2%*%vF
  pval2=pchisq(TF,df=(k-1))
  Alternative2=ifelse(any(vF>1e-06) & any(vF< -1e-06),"Indeterminate","Second Order Stochastic Dominance") 
  Second_Order = list(Statistic = TF,P_value=pval2,Alternative=Alternative2)
  
  #Third-Order Stochastic dominance:
  vC = IF%*%IF%*%If%*%v; vC=matrix(vC,ncol = 1)
  M_inv3 = MASS::ginv(t(IF%*%IF%*%If)%*%moMat%*%(IF%*%IF%*%If))
  TC = t(vC)%*%M_inv3%*%vC
  pval3=pchisq(TC,df=(k-1))
  Alternative3=ifelse(any(vC>1e-06) & any(vC< -1e-06),"Indeterminate","Third Order Stochastic Dominance") 
  Third_Order = list(Statistic=TC,P_value=pval3,Alternative=Alternative3)
  val = list(First_Order=First_Order,Second_Order=Second_Order,Third_Order=Third_Order)
return(val)
}

#----------------------------------------------------
#' @export
stochdomUM <- function(xx) UseMethod("stochdomUM")
#----------------------------------------------------
#' @export
stochdomUM.semipar1<- function(xx){
  X=xx$X; Y=xx$Y
  return(SD_semipar(X,Y)$First_Order)
}
#----------------------------------------------------
#' @export
stochdomUM.semipar2<- function(xx){
  X=xx$X; Y=xx$Y
  return(SD_semipar(X,Y)$Second_Order)
}
#----------------------------------------------------
#' @export
stochdomUM.semipar3<- function(xx){
  X=xx$X; Y=xx$Y
  return(SD_semipar(X,Y)$Third_Order)
}
#----------------------------------------------------
#' @export
stochdomUM.semiparall<- function(xx){
  X=xx$X; Y=xx$Y
  return(SD_semipar(X,Y))
}
#===========================================================================================#
#' @export
cumdist<- function(Y,y){
  ly = length(y)
  val<- rep(NA,ly)
  for(j in 1:ly) {val[j]=length(which(Y<=y[j]))/length(Y)}
  return(val)
}

#----------------------------------------------------
#' @export
stochdomUM.Bootstrap0<- function(xx){
  X=xx$X; Y=xx$Y;B=xx$B
  lX = length(X); lY=length(Y)
  N = min(lX,lY)#shorter vector length
  if(N>=10000){
    k = 512
  }else if(N<36){
    stop("Not enough observations in shorter vector")
  }else{
    k=floor(5.12*sqrt(N))
  }
  Dv = quantile(c(X,Y),probs = seq(0,1,length.out = (k+1))) #k intervals over support of X,Y
  Teq = sqrt(lX*lY/(lX+lY))*max(abs(cumdist(X,Dv)-cumdist(Y,Dv)))
  Z = c(X,Y); lZ = length(Z)

  fn<- function(b){
    set.seed(b); id = sample(1:lZ,lZ,rep=T)
    idx = id[1:lX]; idy = id[(lZ-lY+1):lZ]
    tt=sqrt(lX*lY/(lX+lY))*max(abs(cumdist(Z[idx],Dv)-cumdist(Z[idy],Dv)))
    return(tt)
  }
  Tnb = sapply(1:B,fn); pvalue = length(which(Tnb>Teq))/B
  val = list(statistic=Teq,pvalue=pvalue,alternative="non-equality")
  return(val)
}

#----------------------------------------------------
#' @export
stochdomUM.Bootstrap1<- function(xx){
  X=xx$X; Y=xx$Y;B=xx$B
  lX = length(X); lY=length(Y)
  N = min(lX,lY)#shorter vector length
  if(N>=10000){
    k = 512
  }else if(N<36){
    stop("Not enough observations in shorter vector")
  }else{
    k=floor(5.12*sqrt(N))
  }
  Dv = quantile(c(X,Y),probs = seq(0,1,length.out = (k+1))) #k intervals over support of X,Y
  Tfsd = sqrt(lX*lY/(lX+lY))*max(cumdist(X,Dv)-cumdist(Y,Dv))
  Z = c(X,Y); lZ = length(Z)
  
  fn<- function(b){
    set.seed(b); id = sample(1:lZ,lZ,rep=T)
    idx = id[1:lX]; idy = id[(lZ-lY+1):lZ]
    tt=sqrt(lX*lY/(lX+lY))*max(cumdist(Z[idx],Dv)-cumdist(Z[idy],Dv))
    return(tt)
  }
  Tnb = sapply(1:B,fn); pvalue = length(which(Tnb>Tfsd))/B
  val = list(statistic=Tfsd,pvalue=pvalue,
             Null="First-order stochastic dominance")
  return(val)
}

#----------------------------------------------------
#' @export
stochdomUM.Bootstrap2<- function(xx){
  X=xx$X; Y=xx$Y;B=xx$B
  lX = length(X); lY=length(Y)
  N = min(lX,lY)#shorter vector length
  if(N>=10000){
    k = 512
  }else if(N<36){
    stop("Not enough observations in shorter vector")
  }else{
    k=floor(5.12*sqrt(N))
  }
  Dv = quantile(c(X,Y),probs = seq(0,1,length.out = (k+1))) #k intervals over support of X,Y
  
  hn<- function(j){ pracma::trapz(Dv[1:j],cumdist(X,Dv[1:j])-cumdist(Y,Dv[1:j]))}
  dz = sapply(1:(k+1),hn)
  Tssd = sqrt(lX*lY/(lX+lY))*max(dz)
  Z = c(X,Y); lZ = length(Z)
  
  fn<- function(b){
    set.seed(b); id = sample(1:lZ,lZ,rep=T)
    idx = id[1:lX]; idy = id[(lZ-lY+1):lZ]
    #set.seed(b); idy = sample(1:lY,lY,rep=T);
    hn<- function(j){ pracma::trapz(Dv[1:j],
                                    cumdist(Z[idx],Dv[1:j])-cumdist(Z[idy],Dv[1:j]))}
    dz = sapply(1:(k+1),hn)
    tssd = sqrt(lX*lY/(lX+lY))*max(dz)
    #tt=sqrt(lX*lY/(lX+lY))*max(cumdist(X[idx],Dv)-cumdist(Y[idy],Dv))
    return(tssd)
  }
  Tnb = parsply(1:B,fn); pvalue = length(which(Tnb>Tssd))/B
  val = list(statistic=Tssd,pvalue=pvalue,
             Null="Second-order stochastic dominance")
  return(val)
}



#===========================================================================================#
#' Stochastic Dominance
#'
#' This function tests the stochastic dominance of two variables \code{X} and \code{Y}
#'
#' @param X vector of first variable
#' @param Y vector of second variable
#' @param order for order of test "equality", 1, 2, 3, "all"
#' @param type "Bootstrap" following Abadie (2002) or "semipar" following Anderson (1996)
#' @return val an ndrs x p matrix of coefficient draws, p is the number of OLS coefficients
#' 
#' @examples # Not run
#' # set.seed(30); X=rnorm(200,mean = 0.2,sd=2); set.seed(30); Y=rnorm(400,0.24,sd=1)
#' # stochdom(Y,X,order = 1); stochdom(X,Y,order = 2)
#' # set.seed(30); X=rchisq(1000,df=14); set.seed(30); Y=rnorm(278,mean=6,sd=2)
#' # stochdom(Y,X,order = 1);stochdom(X,Y,order = 1); 
#' # stochdom(X,Y,order = 2); stochdom(Y,X,order = 2)
#' 
#' @export

stochdom<- function(X,Y,order=1,type="Bootstrap",B=2000){
  xx<- classifSD(X,Y,order=order,type=type,B=B)
  return(stochdomUM(xx))
  
}


