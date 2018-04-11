\name{DataRemix}
\alias{DataRemix}
\title{Main DataRemix function}
\usage{
DataRemix(svdres, fn, lower_limit = c(1,-1,0), upper_limit = c(length(transData.svd$d), 1,1), num_of_initialization = 5
                      num_of_thompson = 600, basis = omega, xi = 0.1, full = T, verbose = T,...)
}                      
\arguments{
\item{svdres}{The full SVD decomposition of input gene expression profile X, which is equal to svd(X)}

\item{fn}{The custom function which takes svdres and some additional custom parameters as input and output the biological metric to be maximized}

\item{lower_limit}{The lower bounds for k, p and \eqn{\mu} are set to be 1, -1 and 0 by default}

\item{upper_limit}{The upper bounds for k, p and \eqn{\mu} are set to be the maximum number of loadings in SVD, 1 and 1 by default}

\item{num_of_initialization}{The number of steps of Random Search, which is 5 by default}

\item{num_of_thompson}{The number of steps of Thompson Sampling, which is 600 by default}

\item{basis}{The default mt is set to be 2,500. You can specify a smaller mt by using omega[1:2000,] for example, if you want to set mt to be 2,000.}

\item{xi}{A pre-defined probability which ensures the search doesn't stuck in the local optimum, which is 0.1 by default}

\item{full}{If full is set to be True, the output will be the complete search history. Otherwise the output will only be the optimal combination of k, p and \eqn{\mu} and the corresponding biological metric}

\item{verbose}{If verbose is set to be True, then after each Thompson Sampling iteration, the function will print out the intermediate result.}
}

\value{
\item{para}{If verbose is set to be True, the dimension of para is (num_of_initialization+num_of_thompson)-by-4(k, p, \eqn{\mu}, biological metric). Otherwise, para is only a single vector with length 4}
}
\description{
Main DataRemix function
}
