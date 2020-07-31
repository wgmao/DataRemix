auc_pr<-function(obs, pred) {
  xx.df <- prediction(pred, obs)
  perf  <- performance(xx.df, "prec", "rec")
  xy    <- data.frame(recall=perf@x.values[[1]], precision=perf@y.values[[1]])
  xy <- subset(xy, !is.nan(xy$precision))
  #add a point at 0
  xy <- rbind(c(0, xy[1,2]), xy)
  res   <- trapz(xy$recall, xy$precision)
  return(res)
}#auc_pr

simpleAUC<-function(lab, value){
  value=as.numeric(rank(value)-1)
  posn=as.numeric(sum(lab==1))
  negn=as.numeric(sum(lab!=1))
  stat=sum(value[lab==1])-posn*(posn+1)/2
  return(stat/(posn*negn))
}#simpleAUC

perPathPR<-function(pathPredict, GS){
  length=ncol(pathPredict)
  pathPR=double(length)
  for(i in 1:length){
    x=pathPredict[,i]
    y=GS[,i]
    auc=auc_pr(y, x)
    pathPR[i]=auc
  }
  return(pathPR)
}#perPathPR


perPathAUC<-function(pathPredict, GS){
  length=ncol(pathPredict)
  pathPR=double(length)
  for(i in 1:length){
    x=pathPredict[,i]
    y=GS[,i]
    auc=simpleAUC(y, x)
    pathPR[i]=auc
  }
  return(pathPR)
}#perPathAUC


corMatToAUC=function(data, GS, objective = "mean.AUC"){
  #filter on GS
  ##############################################
  GS.pathway.sum <- apply(GS,2,sum)
  GS <- GS[,which(GS.pathway.sum > 2)]
  
  GS.gene.sum <- apply(GS,1,sum)
  GS <- GS[which(GS.gene.sum > 0), ]
  data <- GS[which(GS.gene.sum > 0), which(GS.gene.sum > 0)]
  
  
  #self-correlation is 0
  ##############################################
  diag(data)=0
  data[is.na(data)]=0
  pathPredict=data%*%GS
  #fix values for genes that are in the pathway since they only get n-1 correlations
  nGenes=colSums(GS)
  #factor to multiply by
  pathwayF=unlist(lapply(nGenes, function(x){x/(x-1)}))
  for(i in 1:ncol(GS)){
    iigenes=which(GS[,i]==1)
    pathPredict[iigenes,i]=  pathPredict[iigenes,i]*pathwayF[i]
  }
  
  PATHAUPR=perPathPR(pathPredict, GS)
  PATHAUC=perPathAUC(pathPredict, GS)
  
  if (objective == "mean.AUC"){
    return(c(mean(PATHAUPR),mean(PATHAUC)))    
  }else if (objective == "mean.AUPR"){
    return(c(mean(PATHAUC), mean(PATHAUPR)))
  }else if (objective == "median.AUC"){
    return(c(median(PATHAUPR),median(PATHAUC)))    
  }else{
    return(c(median(PATHAUC), median(PATHAUPR)))
  }#else

}#corMatToAUC



HCP=function(data, covariates, k, L1, L2, L3, max.iter, trace=F, return.all=F){
  
  Y=scale(t(data))/sqrt(ncol(data)-1)
  C=scale(covariates)/sqrt(nrow(covariates)-1)
  
  
  nc=ncol(covariates)
  ng=nrow(data)
  ns=ncol(data)
  
  
  U=matrix(0,nrow=nc, ncol=k )
  Z=matrix(0, nrow=ns, ncol=k)
  svdres=svd(t(Y))
  B=t(svdres$u[, 1:k])
  
  for ( i in 1:max.iter){
    err0=sum((Y-Z%*%B)^2)+sum((Z-C%*%U)^2)*L1+sum(B^2)*L2+sum(U^2)*L3
    if(trace){
      print(err0)
    }
    Z=(Y%*%t(B)+L1*C%*%U)%*%solve(B%*%t(B)+L1*diag(k))
    B=solve(t(Z)%*%Z+L2*diag(k))%*%t(Z)%*%Y
    U=solve(t(C)%*%C*L1+L3*diag(nc))%*%t(C)%*%Z*L1    
  }
  
  if(return.all){
    return(list(res=(Y-Z%*%B), B=B, Z=Z, U=U))
  }
  else{
    return(t(Y-Z%*%B))
    
  }
}
