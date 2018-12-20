#library(MASS)


tr <- function(matrix){
  return(sum(diag(matrix)))
}#tr


normF <- function(x){
  return(sqrt(sum(x^2)))
}#fastcor


kernel <- function(x,y){
  #x <- as.matrix(x)
  #y <- as.matrix(y)
  #if (ncol(x)==1){
  #  x <- t(x)
  #}#if
  #if (ncol(y)==1){
  #  y<-t(y)
  #}#if
  res <- matrix(0, nrow = nrow(y), ncol = nrow(x))
  for (i in 1:nrow(y)){
    for (j in 1:nrow(x)){
      res[i,j] <- exp(-0.5*normF(y[i,]-x[j,]))
    }#for j
  }#for i
  return(res)
}#kernel



marginal_likelihood <- function(sigma, K, y){
  delta_matrix <- diag(sigma^2, nrow = nrow(K))
  ml <- -0.5*log(det(K+delta_matrix))-0.5*t(y)%*%solve(K+delta_matrix)%*%y
  return(-ml)
}#marginal likelihood



feature_map <- function(x, basis, b){
  x_row <- nrow(x) #num of points
  x_col <- ncol(x) #dimension
  mt <- nrow(basis)
  res <- matrix(1, nrow = mt, ncol = x_row) #mt: number of bases
  
  
  for (j in 1:x_row){
    for (i in 1:nrow(res)){
      res[i,j] <- sqrt(2/mt)*cos(sum(x[j,]*basis[i,])+b[i])
    }#for i
  }#for j
  
  return(res)
}#feature_map



bayesian <- function(x, theta, basis, b){
  dim_x <- length(x)
  x <- matrix(x, nrow = 1)
  phi <- feature_map(x, basis, b)
  
  approx <- t(phi)%*%theta
  return(-approx)
}#bayesian



bayesian_A <- function(record, y, sigma, basis, b){
  phi_rec <- feature_map(record, basis, b)
  
  #posterior infer
  A <- 1/sigma^2*phi_rec%*%t(phi_rec)+diag(1, nrow = nrow(basis))
  A_inv <- solve(A)
  mu <- 1/sigma^2*A_inv%*%phi_rec%*%y
  variance <- A_inv
  #sample from posterior
  theta <- mvrnorm(1, mu = mu, Sigma = variance)
  theta <- matrix(theta,ncol = 1)
  
  return(theta)
}#bayesian



thompson_sampling <- function (lower, upper, record, K, y, basis, b) 
{
  #infer sigma
  set.seed(1)
  #sigma_res <- optim(0,marginal_likelihood,method = "Brent", lower=-10, upper = 10, K=K,y=y)
  #sigma <- sigma_res$par
  sigma <- 0.1
  
  #infer the bayesian
  theta <- bayesian_A(record, y, sigma, basis, b)
  
  par <- c()
  minimum <- 0
  seed_num <- 20
  sp_list1 <- runif(seed_num, min = lower[1], max = upper[1])
  sp_list2 <- runif(seed_num, min = lower[2], max = upper[2])
  sp_list3 <- runif(seed_num, min = lower[3], max = upper[3])
  sp_list <- cbind(sp_list1,sp_list2,sp_list3)
  for (i in 1:nrow(sp_list)){
    sp <- sp_list[i,]
    bayesian_res <- optim(sp, bayesian, method = "L-BFGS-B", lower = lower, upper= upper, theta = theta, basis = basis, b = b)
    if (bayesian_res$value < minimum){
      minimum <- bayesian_res$value
      par <- bayesian_res$par
    }#if
  }#for i
  
  return(list(par=par))
}#thompson_sampling




main <- function(record, lower_limit, upper_limit, basis, b){
  #log scale on mu
  record[,length(lower_limit)] <- log10(record[,length(lower_limit)])
  lower_limit[length(lower_limit)]<- log10(lower_limit[length(lower_limit)])
  upper_limit[length(upper_limit)] <- log10(upper_limit[length(upper_limit)])
  
  #scale
  mean_val <- apply(record,2,mean)
  sd_val <- apply(record,2,sd)
  for (j in 1:length(sd_val)){
    if (sd_val[j]==0){
      sd_val[j] <-1
      mean_val[j] <- 0
    }#if
  }#for j
  record_scale <- sweep(record,2,mean_val,"-")
  record_scale <- sweep(record_scale,2,sd_val,"/")
  
  
  #pre-computation
  y <- matrix(record_scale[,ncol(record_scale)],ncol=1)
  record <- record_scale[,-ncol(record_scale)]
  K <- kernel(record, record)
  
  
  #scale the limit
  lower_scale <- c()
  upper_scale <- c()
  for (j in 1:length(lower_limit)){
    lower_scale <- c(lower_scale, (lower_limit[j]-mean_val[j])/sd_val[j] )
    upper_scale <- c(upper_scale, (upper_limit[j]-mean_val[j])/sd_val[j] )
  }
  
  res <- thompson_sampling(lower = lower_scale, upper = upper_scale, record, K, y, basis, b) 
  
  
  #evaluate the new x
  new_comb <- c()
  for (j in 1:length(lower_limit)){
    if (j==1){
      new_comb <- c(new_comb, min(max(round(res$par[j]*sd_val[j]+mean_val[j]),lower_limit[j]), upper_limit[j])  )  
    }else if (j==3){
      new_comb <- c(new_comb, 10^min(max(res$par[j]*sd_val[j]+mean_val[j],lower_limit[j]), upper_limit[j])   )  
    }else{
      new_comb <- c(new_comb, min(max(res$par[j]*sd_val[j]+mean_val[j],lower_limit[j]), upper_limit[j])   )  
    }#else
  }#for j
  return(new_comb)
}#main



random_sample <- function(lower, upper){
  para1 <- sample(lower[1]:upper[1],1)
  para2 <- runif(1, min = lower[2], max = upper[2])
  para3 <- runif(1, min = lower[3], max = upper[3])
  return(c(para1,para2,para3))
}#random_sample


coarsed_grid <- function(lower, upper){
  loc <- c(0.25, 0.5, 0.75)
  para <- c()
  for (i in 1:length(loc)){
    for (j in 1:length(loc)){
      for (k in 1:length(loc)){
         para1 <- floor(lower[1]+loc[i]*(upper[1]-lower[1]))
         para2 <- lower[2]+loc[j]*(upper[2]-lower[2])
         para3 <- 10^(log10(lower[3])+loc[k]*(log10(upper[3])-log10(lower[3])))
         para <- rbind(para, c(para1, para2, para3))
      }#for k
    }#for j
  }#for i
  return(para)
}#coarsed_grid



SVDcombine<-function(svdData, k, power, lambda=0){
  if (k >1){
    SVDcombine=svdData$u[, 1:k] %*% diag(svdData$d[1:k]^(power)) %*% t(svdData$v[, 1:k])  
  }else{
    SVDcombine= (svdData$d[1:k]^(power)) *matrix(svdData$u[, 1:k],ncol = 1) %*% matrix(svdData$v[, 1:k], nrow = 1)  
  }
  
  if(lambda>0){
    nc=ncol(svdData$u)
    if (k ==nc-1){
      SVDcombine=SVDcombine+svdData$d[(k+1):nc]* lambda* matrix(svdData$u[, (k+1):nc],ncol=1) %*% matrix(t(svdData$v[, (k+1):nc]),nrow=1)
    }else if (k < nc-1){
      SVDcombine=SVDcombine+lambda*svdData$u[, (k+1):nc] %*% diag(svdData$d[(k+1):nc]) %*% t(svdData$v[, (k+1):nc])  
    }#else if
  }#if
  rownames(SVDcombine)=rownames(svdData$u)
  colnames(SVDcombine)=colnames(svdData$u)
  SVDcombine
}#SVDcombine


#DESCRIPTION: LazyData: true
DataRemix <- function(svdres, fn, k_limits = c(1, ceiling(length(svdres$d)/2)), p_limits = c(-1,1), mu_limits = c(1e-12,1), num_of_initialization = 5, num_of_thompson = 600, basis = omega, xi = 0.1, full = T, verbose = T, ...){
    mt <- nrow(basis)
    set.seed(1)
    b <- runif(mt)*2*pi

    lower_limit <- c(k_limits[1], p_limits[1], mu_limits[1])
    upper_limit <- c(k_limits[2], p_limits[2], mu_limits[2])
    
    #Initialization
    history <- c()
    if (num_of_initialization==0){
        para <- coarsed_grid(lower_limit, upper_limit)
        for (i in 1:nrow(para)){
          rec.term <- fn(SVDcombine(svdres, k= para[i,1], p=para[i,2], lambda = para[i,3]), ...)
          history <- rbind(history,c(para[i,], rec.term))
        }#for i
    }else{
    for (i in 1:num_of_initialization){
      para <- random_sample(lower_limit, upper_limit)
      rec.term <- fn(SVDcombine(svdres, k= para[1], p=para[2], lambda = para[3]), ...)
      history <- rbind(history,c(para, rec.term))
    }#for i
    }#Random Search
    
    
    #Thompson
    record <- history[,c(1:length(lower_limit),ncol(history))]
    
    for (i in 1:num_of_thompson){
      uniform <- runif(1)
      if (uniform < xi){
        para <- random_sample(lower_limit, upper_limit)
      }else{
        para <- main(record,lower_limit, upper_limit, basis, b) 
      }#else
      
      #duplicated para
      while (length(which(record[,1]==para[1] & record[,2]==para[2] & record[,3] == para[3])) > 0){
        para <- random_sample(lower_limit, upper_limit)
      }#while
      
      rec.term <- fn(SVDcombine(svdres, k= para[1], p=para[2], lambda = para[3]), ...)
      history <- rbind(history,c(para, rec.term))
      record <- rbind(record, c(para, rec.term[length(rec.term)]))
      if (verbose){
        print(c(i, para, rec.term[length(rec.term)]))  
      }#if
    }#for i
    
    #output record or only the best
    if (full){
      return(list(para = record, full = history))
    }else{
      index = which.max(record[ncol(record)])
      return(list(para = record[index,], full = history[index,]))
    }#else
}#DataRemix


