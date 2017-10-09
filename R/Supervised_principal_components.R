#===========================================================================================#
#' @export

sigf<- function(j,y,x,lev1){ # is x significant at level lev of significance? Returns logical
  cc<-summary(stats::lm(y~x[,j]))$coefficients[-1,4] <lev1
  return(cc)
}
#===========================================================================================#
#' Significant feature selection
#'
#' \code{sfeat} produces a subset a matrix x with features significant at at most \code{lev1}
#' level of significance.
#'
#' @param y predictand -  a vector of length n
#' @param x n x p matrix, with possibly p>n
#' @param lev1 maximum level of significance to consider; (default 0.1)
#' @param minft minimum number of features to select. If the number of selected features is
#' less than \code{minft}, return the entire design matrix \code{x}; (default 6)
#' @return xslct a subset matrix of \code{x} whose columns comprises selected features
#' @return slctfeat a vector of indices of selected features
#'
#' @export

sfeat<- function(y,x,lev1=0.1,minft=6){ # select features for principal components
  vc<- 1:ncol(x);xs<- scale(x)
  cc<- sapply(vc,sigf,y=y,x=xs,lev1=lev1)
  slctfeat<- which(cc)
  nslct<- length(slctfeat)
  if(nslct<minft){
    slctfeat<- vc
  }
  xslct<- as.matrix(x[,slctfeat])
  return(list(xslct=xslct,slctfeat=slctfeat))
}
#===========================================================================================#

# Example: dd<- sfeat(y,x)


#===========================================================================================#
#' Supervised principal component
#'
#' This function conducts principal components on a subset matrix \code{xx}
#'
#' @param dat data set with first column being the predictand/dependent variable
#' @param xvec predictor vector of interest. Its supervised component form will be the first row of the matrix returned
#' @param lev1 lev1 maximum level of significance to consider; (default 0.1)
#' @param lev2 highest fraction of cumulative variation of supervised principal components: default 0.95
#' @param maxc maximum number of supervised principal components to return. Note, this option overrides
#' the \code{lev} option if \code{maxc} is set to a number smaller than the one dictated by \code{lev2}
#' @param minft minimum number of features to select. If the number of selected features is
#' less than \code{minft}, return the entire design matrix \code{x}; (default 6)
#' @return slctcomp  a matrix of maxc supervised principal components.
#'
#' @examples
#' dat<- prdtrade::DAT_ALL[1:75,]; dat<- Ddat[-20,]; xvec=Ddat[20,-1];
#' spc<- supc(dat = dat,xvec=xvec)
#'
supc<- function(dat,xvec,lev1=0.1,lev2=0.95,maxc=3, minft=6){
  y=dat[,1]; x=dat[,-1]
  dd<-sfeat(y,x,lev1 = lev1,minft = minft)
  xx<-dd$xslct; slctfeat<- dd$slctfeat

  if(!is.null(xvec)){
  xx<- as.matrix(rbind(xvec[slctfeat],xx))
  }
  rpca <- stats::prcomp(xx, scale = TRUE); dc<-(rpca$sdev/sum(rpca$sdev))
  csdc<- cumsum(dc); th<-min(which(csdc>lev2))
  if(th>maxc){
    th=maxc
  print.noquote(paste("Selecting only maxc=",maxc," principal components."))
  }
  slctcomp<- as.matrix(rpca$x[,1:th])
  return(slctcomp)
}





