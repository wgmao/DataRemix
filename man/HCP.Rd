\name{HCP}
\alias{HCP}
\title{HCP}
\usage{
  HCP(data, covariates, k, L1, L2, L3, max.iter, trace=F, return.all=F)
}                      
\arguments{
  \item{data}{The data to be processed}
  
  \item{covariates}{Covariate matrix associated with samples}
  
  \item{k}{Specify number of hidden factors}
  
  \item{L1}{Model parameter-penalty on \eqn{||Z-CU||_2^2}. Here \eqn{C} denotes the scaled covariates}
  
  \item{L2}{Model parameter-penalty on \eqn{||B||_2^2}}
  
  \item{L3}{Model parameter-penalty on  \eqn{||U||_2^2}}
  
  \item{max.iter}{Maximum number of iterations to perform}
  
  \item{trace}{Displace progress information or not}
  
  \item{return.all}{Return all reuslts or only the residual}
  
}

\value{
\item{res}{Residual out of decomposition}
\item{B}{Latent factors (if return.all=T)}
\item{Z}{Loading matrix (if return.all=T)}
\item{U}{Coefficient matrix for covariates (if return.all=T)}
}
\description{
Hidden Covariates with Prior (HCP) implementation
}
