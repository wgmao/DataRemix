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

corMatToAUC=function(data, GS){
  #self-correlation is 0
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
  return(c(mean(PATHAUPR),mean(PATHAUC)))
}#corMatToAUC
