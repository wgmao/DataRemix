library(DataRemix)

#define fn
reconstruct <- function(X_reconstruct, X, penalty){
  return(-sum((X-X_reconstruct)^2)+penalty)
}#reconstruct


#generate X and svdres
set.seed(1)
num_of_row <- 100
num_of_col <- 9
X <- matrix(rnorm(num_of_row*num_of_col), nrow = num_of_row, ncol = num_of_col)
svdres <- svd(X)

#load basis
basis_short <- basis[1:2000,]


#start the optimization
DataRemix.res <- DataRemix(svdres, reconstruct, lower_limit = c(1,-1,0), upper_limit = c(length(svdres$d), 1,1), num_of_initialization = 5,                   num_of_thompson = 50, basis = basis_short, xi = 0.1, full = T, verbose = T, X = X, penalty = 100)