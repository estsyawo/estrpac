% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Distri_Sim_functions.R
\name{Band.Fn.objs}
\alias{Band.Fn.objs}
\title{\code{Band.Fn.objs} constructs Chernozhukov et al. 2013 Uniform Confidence Bands}
\usage{
Band.Fn.objs(Fn, Boot.Fn, n)
}
\arguments{
\item{Fn}{G-length vector of the function of interest}

\item{Boot.Fn}{a \eqn{G\times L} matrix of samples of the function Fn}

\item{n}{sample size of the data used to compute Fn}
}
\value{
Sigv:

tbv:
}
\description{
\code{Band.Fn.objs} constructs Chernozhukov et al. 2013 Uniform Confidence Bands
}
\examples{
set.seed(1); Band.Fn.objs(Fn=rep(0,5),Boot.Fn=matrix(rnorm(1000),nrow=5),n=100)
}
