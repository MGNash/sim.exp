getcdf = function(nsubj,nrep){
  corlist = rep(NA,nrep)
  for(i in 1:nrep){
    x = 1:nsubj
    y = sample(x,nsubj,replace=FALSE)
    corlist[i] = cor(x,y,method="spearman")
  }
  return(ecdf(corlist))
}

sim.exp=function(nrep=1,nsubj = 16,nprot=600,ntrue=30,FDR=.5,rho=.5,spearman=TRUE){
  cdf.exp=getcdf(nsubj,10000)
  tpos = rep(NA,nrep)
  fpos = rep(NA,nrep)
  pval.i = function(y,x,spearman=TRUE){
    #x = xy[,1]
    #y = xy[,2]
    if(spearman){
      corx = cor(x,y,method="spearman")
      pval = 2*min(cdf.exp(corx),(1-cdf.exp(corx)))
      return(pval)
    }
    else {
      #return(pvalue(spearman_test(y~x)))
      return(cor.test(x,y)$p.value)
    }
  }
  for(i in 1:nrep){
    x = rnorm(nsubj)
    y.1 = matrix(rep(x,each=ntrue)+rnorm(ntrue*nsubj,0,sqrt((rho^-2)-1)),nrow=ntrue)
    y.2 = matrix(rnorm((nprot-ntrue)*nsubj),ncol=nsubj)
    y = rbind(y.1,y.2)
    pval.xy = apply(y,1,pval.i,x=x,spearman=spearman)
    jitter = rep(0,nprot)
    if(spearman){jitter = runif(nprot,0,.00001)}
    pval.xy.j = pval.xy + jitter
    ptable = data.frame(pval = pval.xy.j, ID = 1:nprot, rank = NA)
    ptable = ptable[order(ptable$pval,decreasing=FALSE),]
    ptable$rank=1:nprot
    BMcv = (ptable$rank/nprot)*FDR
    underCV = ptable$pval < BMcv
    if(sum(underCV)>=1){
      max.sig = which.max(ptable$rank[underCV])
      ptable$sig = ptable$rank <= max.sig
      tpos[i] = sum(ptable$sig & (ptable$ID <= ntrue))
      fpos[i] = sum(ptable$sig & (ptable$ID > ntrue))
    } else {
      tpos[i] = 0
      fpos[i] = 0
    }
  }
  return(data.frame(true_positive = tpos,false_positive = fpos))
}