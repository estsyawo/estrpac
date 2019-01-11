#===========================================================================================#
#' Model Selection
#'
#' This function generates a best subset selection matrix of zeros and ones corresponding to
#' inclusion/exclusion of covariates in the data frame \code{daten}.
#'
#' @param daten the data frame of covariates.
#' @param sbset if specified as a number from 1 to number of columns of \code{daten}, the matrix
#' returned corresponds to selections with at least \code{sbset} covariates. If a negative number
#' from -1 to negative the number of columns, the matrix returned corresponds to selections with
#' at most \code{sbset} covariates.
#' @return gridm the matrix grid whose rows correspond to inclusion (1) and exclusion (0) of
#' covariates.
#'
#' @examples
#' gz<- matrix(rnorm(20),ncol = 4); grd<- modslct(gz)
#' gz[,which(grd[8,]==1)] ## Use 8th subset
#' (grd1<- modslct(gz,3));(grd1<- modslct(gz,-2))
#'
#' @export
modslct<- function(daten,sbset=NULL){
  nc = ncol(daten)
  gz<-as.matrix(rbind(rep(0,nc),rep(1,nc)))
  lgz<- split(gz, rep(1:ncol(gz), each = nrow(gz)))
  gridm<-expand.grid(lgz)

  if(!is.null(sbset)){
    suma<- apply(gridm,1,sum)
    gridm<- gridm[sign(sbset)*suma>=sbset,]
  }
  return(gridm)
}
#===========================================================================================#
#' Probability value of Sharpe Ratio
#'
#' This function computes the probaility value of the Sharpe ratio
#'
#' @param ydist a vector of draws from the distribution of excess returns
#' @param yt the value of comparison, i.e. we test if expected excess returns is equal to yt; (default; yt=0)
#' @param yhat expected excess returns; provide this option if you do not have the distribution \code{ydist}
#' @param ysd the volatility of excess returns; provide this input if you do not have the distribution \code{ydist}
#' @return pval the probability value
#'
#' @examples
#' d<- rnorm(5000,-0.023,1.3); testSR(d); testSR(yt=0.5,yhat=0.11,ysd=1.01)
#'
#' @export

testSR<- function(ydist=NULL,yt=0,yhat=NULL,ysd=NULL){
  if(!is.null(ydist)){
    yhat<- mean(ydist); ysd = sd(ydist)
  }
  p<- (yhat-yt)/ysd
  SR<- yhat/ysd
  sdSR<- sqrt(1+0.5*(SR^2))
  z<- p/sdSR
  pval<-2*pnorm(-abs(z)) #fail to reject null hypothesis
  return(pval)
}
#===========================================================================================#
#' Parallel compute
#' 
#' This function takes a function and margin and parallel computes, returning output in a simplified form.
#' 
#' @param X a vector as used in lapply() family of functions
#' @param fn the function to be run in parallel
#' @param type the type of parallel computation; \code{"PSOCK"} which is compatible with all operating
#' systems or \code{"FORK"} which is only compatible on Mac/Linux platforms
#' #' @param nc number of clusters to use
#' @param ... additional inputs for the function \code{fn}.
#' @return val output in a simplified form: vector or matrix
#' 
#' @export
parsply<- function(X,fn,type="PSOCK",nc=1,...){#parallelise a function
  c1<-parallel::makeCluster(nc,type = type)
  val<-parallel::parSapply(c1,X,fn,...)
  parallel::stopCluster(c1) 
  return(val)
}

#===========================================================================================#
#' Parallel compute
#' 
#' This function takes a function and margin and parallel computes, returning output in a list.
#' 
#' @param X a vector as used in lapply() family of functions
#' @param fn the function to be run in parallel
#' @param type the type of parallel computation; \code{"PSOCK"} which is compatible with all operating
#' systems or \code{"FORK"} which is only compatible on Mac/Linux platforms
#' @param nc number of clusters to use
#' @param ... additional inputs for the function \code{fn}.
#' @return val output in a simplified form: vector or matrix
#' 
#' @export
parlply<- function(X,fn,type="PSOCK",nc=1,...){#parallelise a function
  c1<-parallel::makeCluster(nc,type = type)
  val<-parallel::parLapply(c1,X,fn,...)
  parallel::stopCluster(c1) 
  return(val)
}
