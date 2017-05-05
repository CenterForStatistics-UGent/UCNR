### this file contains auxiliary functions needed for reproduction of the simulation studies

unified.design.matrix<-function(miRNA,subject,group,nrMi,nrS,intercept=TRUE,sub.sum=F,UV=NULL) {
  contrasts(miRNA)<-contr.sum
  if (sub.sum) {
    contrasts(subject)<-contr.sum;
  }
  dim.tot<-2*nrMi+nrS-1
  m<-model.matrix(~miRNA+subject+miRNA:group)
  if (!is.null(UV)) {
    X<-(1-UV)
    m<-model.matrix(~miRNA+subject:X+miRNA:group)
    dim.tot<-2*nrMi+nrS
  }
  mm<-m
  for (i in 1:(nrMi-1)) {
    mm[,dim.tot-i]<- mm[,dim.tot-i]-mm[,dim.tot]
  }
  mm<-mm[,-dim.tot]
  if (!intercept) mm<- mm[,-1];
  return(mm);
}

mestdagh.design.matrix<-function(miRNA,group,nrMi,nrS,intercept=TRUE,miRNA.sum=T) {
  if (miRNA.sum) {
    contrasts(miRNA)<-contr.sum
  }
  m<-model.matrix(~miRNA+miRNA:group)
  mmm<-m
  dim.tot<-2*nrMi+nrS-1
  dim.tot.subj<-dim.tot-(nrS-1)
  for (i in 1:(nrMi-1)) {
    mmm[,dim.tot.subj-i]<- mmm[,dim.tot.subj-i]-mmm[,dim.tot.subj]
  }
  mmm<-mmm[,-dim.tot.subj]
  if (!intercept) mmm<- mmm[,-1];
  return (mmm);
}

normalise.qpcr.data<-function (qpcr.data.matrix,LOD=40,nMi,nEach) {
  sample.norm<-matrix(0,nc=ncol(qpcr.data.matrix),nr=nrow(qpcr.data.matrix))
  sample.norm.max<-matrix(0,nc=ncol(qpcr.data.matrix),nr=nrow(qpcr.data.matrix))
  sample.norm.mu<-matrix(0,nc=ncol(qpcr.data.matrix),nr=nrow(qpcr.data.matrix))
  
  for (s in 1:ncol(qpcr.data.matrix)) {
    norm.ct<-c()
    beta<-c()
    mu<-mean(qpcr.data.matrix[,s]);
    for (j in 1:(2*nEach)) {
      expr.x<-c()
      teller<-1;
      for (i in 1:nMi) {
        cq<-as.double(qpcr.data.matrix[(i-1)*(2*nEach)+j,s])
        if (cq!=LOD) {		
          expr.x[teller]<-cq;
          teller<-teller+1;		
        }
      }
      norm.ct[j]<-mean(expr.x);
      beta[j]<-(norm.ct[j]-mu)*(teller/nMi)
    }
    #Normalize
    for (j in 1:(2*nEach)) {
      for (i in 1:nMi) {
        sample.norm[(i-1)*(2*nEach)+j,s]<-as.double(qpcr.data.matrix[(i-1)*(2*nEach)+j,s])-norm.ct[j]
      }
    }
    sample.norm.max[,s]<-sample.norm[,s]
    sample.norm.mu[,s]<-sample.norm[,s]
    
    for (i in 1:nMi) {
      max.mi<-max(as.double(sample.norm[((i-1)*(2*nEach)+1):((i)*(2*nEach)),s]))
      #max.mi<-max(as.double(sample.norm[((i-1)*(2*nEach)+idx.expr),s]))
      for (j in 1:(2*nEach)) {
        cq<-as.double(qpcr.data.matrix[(i-1)*(2*nEach)+j,s])
        if (cq==LOD) {
          sample.norm.max[(i-1)*(2*nEach)+j,s]<-max.mi;
          sample.norm.mu[(i-1)*(2*nEach)+j,s]<-LOD-mu;
        }
        else {
          sample.norm.mu[(i-1)*(2*nEach)+j,s]<-as.double(qpcr.data.matrix[(i-1)*(2*nEach)+j,s])-mu-beta[j];									
        }
      }
    }
    
  }
  return(list(qpcr.norm=sample.norm,qpcr.norm.max=sample.norm.max,qpcr.norm.mu=sample.norm.mu,norm.ct=norm.ct));
}


plot.p<-function(miR,minY=-3,maxY=3,idx=ncol(de.tobit.est),minP=NULL,maxP=NULL) {
  nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow=TRUE), respect=TRUE)
  cols<-c("green","red","blue")#colors()[sample(100:150,5)]"gray","wheat",
  pchs<-c(21,22,23)#sample(19:25,3)  
  rect.length=2

  plot(de.tobit.est.varL[miR,1:idx],col=cols[1],ylim=c(minY,maxY),main=paste("miRNA",miR," (final # UV : ",sum(where.UV[miR,]==1),"(MNA) /",sum(where.UV[miR,]==2)," (MNSC))"),ylab=expression(paste("Estimates Differential Expression ",hat(delta)[i])),xlab="#Censored obs in dataset",type="l",lwd=3,lty=1)
  lines(de.mestdagh.est.varL[miR,1:idx],col=cols[2],lwd=3,lty=2)
  lines(de.mestdagh.est.max[miR,1:idx],col=cols[3],lwd=3,lty=3)
  
  abline(h=de.tobit.orig[miR])
  legend("topleft",c("UCNR","LOD","MNV"),col=cols,lty=1:3,lwd=3)
  grid(NA,4,lwd=2)
  tcks.1<-which(where.UV[miR,]==1)
  tcks.2<-which(where.UV[miR,]==2)
  if (length(tcks.1>0)) {
    #rect(tcks.1-rect.length,minY,tcks.1+rect.length,maxY,lty=1)
    #text(tcks.1)
    points(tcks.1,rep(minY,length(tcks.1)),pch=21,bg="black")
  }
  if (length(tcks.2>0)) {
    #		rect(tcks.2-rect.length,minY,tcks.2+rect.length,maxY,lty=2)
    points(tcks.2,rep(minY,length(tcks.2)),pch=24,bg="grey")
    
  }
  minY<-0
  maxY<-1
  all.ps<-c(p.tobit.est.varL[miR,1:idx],p.mestdagh.est.varL[miR,1:idx],p.mestdagh.est.max[miR,1:idx])
  if (is.null(minP)) {
    minP<-min(all.ps)-0.01
  }
  if (is.null(maxP)) {
    maxP<-1
  }
  plot(p.tobit.est.varL[miR,1:idx],type="l",col=cols[1],ylim=c(minP,maxP),main=paste("miRNA",miR," (# UV : ",sum(where.UV[miR,]==1),"/",sum(where.UV[miR,]==2),")"),ylab="P-values Differential Expression Hypothesis",xlab="#Censored obs in dataset",lty=1,lwd=3)
  lines(p.mestdagh.est.varL[miR,1:idx],type="l",col=cols[2],lty=2,lwd=3)
  lines(p.mestdagh.est.max[miR,1:idx],type="l",col=cols[3],lty=3,lwd=3)
  
  des<-c(de.tobit.est.varL[miR,1:idx],de.mestdagh.est.varL[miR,1:idx],de.mestdagh.est.max[miR,1:idx])
  nfs<-gl(3,idx,labels=c("UCNR","LOD","MNV"))
  boxplot(des~nfs,main=expression(paste("Box plot of ",hat(delta)[i])))
  abline(h=de.tobit.orig[miR])
}


plot.p.extra<-function(miR,minY=-3,maxY=3,idx=ncol(de.tobit.est),minP=NULL,maxP=NULL) {
  nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow=TRUE), respect=TRUE)
  cols<-c("green","red","blue","black","grey","yellow")#colors()[sample(100:150,5)]"gray","wheat",
  pchs<-c(21,22,23,24,25,26)#sample(19:25,3)  
  rect.length=2
  plot(de.tobit.est.varL[miR,1:idx],col=cols[1],ylim=c(minY,maxY),main=paste("miRNA",miR," (final # UV : ",sum(where.UV[miR,]==1),"(MNA) /",sum(where.UV[miR,]==2)," (MNSC))"),ylab=expression(paste("Estimates Differential Expression ",hat(delta)[i])),xlab="#Censored obs in dataset",type="l",lwd=3,lty=1)
  lines(de.mestdagh.est.varL[miR,1:idx],col=cols[2],lwd=3,lty=2)
  lines(de.mestdagh.est.max[miR,1:idx],col=cols[3],lwd=3,lty=3)
  lines(de.mestdagh.est.knn[miR,1:idx],col=cols[4],lwd=3,lty=4)
  lines(de.tobit.est.weighted[miR,1:idx],col=cols[5],lwd=3,lty=5)
  lines(de.mestdagh.mod.glob[miR,1:idx],col=cols[6],lwd=3,lty=6)
  
  abline(h=de.tobit.orig[miR])
  legend("topleft",c("UCNR","LOD","MNV","KNN","UCNRW","GLOBMOD"),col=cols,lty=1:6,lwd=3)
  grid(NA,4,lwd=2)
  tcks.1<-which(where.UV[miR,]==1)
  tcks.2<-which(where.UV[miR,]==2)
  if (length(tcks.1>0)) {
    points(tcks.1,rep(minY,length(tcks.1)),pch=21,bg="black")
  }
  if (length(tcks.2>0)) {
    points(tcks.2,rep(minY,length(tcks.2)),pch=24,bg="grey")
  }
  minY<-0
  maxY<-1
  all.ps<-c(p.tobit.est.varL[miR,1:idx],p.mestdagh.est.varL[miR,1:idx],p.mestdagh.est.max[miR,1:idx],p.mestdagh.est.knn[miR,1:idx],p.tobit.est.weighted[miR,1:idx],p.mestdagh.mod.glob[miR,1:idx])
  if (is.null(minP)) {
    minP<-min(all.ps)-0.01
  }
  if (is.null(maxP)) {
    maxP<-1
  }
  plot(p.tobit.est.varL[miR,1:idx],type="l",col=cols[1],ylim=c(minP,maxP),main=paste("miRNA",miR," (# UV : ",sum(where.UV[miR,]==1),"/",sum(where.UV[miR,]==2),")"),ylab="P-values Differential Expression Hypothesis",xlab="#Censored obs in dataset",lty=1,lwd=3)
  lines(p.mestdagh.est.varL[miR,1:idx],type="l",col=cols[2],lty=2,lwd=3)
  lines(p.mestdagh.est.max[miR,1:idx],type="l",col=cols[3],lty=3,lwd=3)
  lines(p.mestdagh.est.knn[miR,1:idx],type="l",col=cols[4],lty=4,lwd=3)
  lines(p.tobit.est.weighted[miR,1:idx],type="l",col=cols[5],lty=5,lwd=3)
  lines(p.mestdagh.mod.glob[miR,1:idx],type="l",col=cols[6],lty=6,lwd=3)
  
  des<-c(de.tobit.est.varL[miR,1:idx],de.mestdagh.est.varL[miR,1:idx],de.mestdagh.est.max[miR,1:idx],de.mestdagh.est.knn[miR,1:idx],de.tobit.est.weighted[miR,1:idx],de.mestdagh.mod.glob[miR,1:idx])
  nfs<-gl(6,idx,labels=c("UCNR","LOD","MNV","KNN","UCNRW","MODGLOB"))
  boxplot(des~nfs,main=expression(paste("Box plot of ",hat(delta)[i])))
  abline(h=de.tobit.orig[miR])
}

plot.nc<-function(subj,minY=NULL,maxY=NULL,idx=ncol(tobit.beta.est)) {
  par(mfrow=c(1,1))
  miRNA<-factor(testData$miRNA)
  nMi<-length(unique(miRNA))
  cols<-c("green","red","blue")#colors()[sample(100:150,5)]
  pchs<-c(21,22,23)#sample(19:25,3)  
  if (is.null(minY)) {
    minY<-min(c(tobit.beta.est.varL[subj,1:idx],norm.est.ct[subj,1:idx]),na.rm=T)
  }
  if (is.null(maxY)) {
    maxY<-max(c(tobit.beta.est.varL[subj,1:idx],norm.est.ct[subj,1:idx]),na.rm=T)
  }
  plot(tobit.beta.est.varL[subj,1:idx],col=cols[1],pch=pchs[1],ylim=c(minY,maxY),main=paste("Sample ",subj),ylab="Estimates Normalization Constant",xlab="#Censored obs in dataset",type="l",lty=1,lwd=3)
  lines(norm.est.ct[subj,1:idx],col=cols[2],pch=pchs[2],lty=2,lwd=3)
  abline(h=tobit.beta[subj])
  legend("topleft",c("Latent Mean Normalization","Global Mean Normalization"),col=cols,lwd=3,lt=1:2)
}

plot.mse.mir<-function(miR,minY=0,maxY=4,idx=ncol(tobit.beta.est)) {
  par(mfrow=c(1,1))
  cols<-c("green","red","blue")#colors()[sample(100:150,5)]
  pchs<-c(21,22,23)#sample(19:25,3)  
  tcks.1<-which(where.UV[miR,]==1)
  tcks.2<-which(where.UV[miR,]==2)
  plot(mse.de.tobit.varL[miR,1:idx],col=cols[1],pch=pchs[1],type="l",ylim=c(minY,maxY),main="MSE per amount of censored obs",ylab="MSE",xlab="#Censored obs in dataset")
  lines(mse.de.mestdagh.varL[miR,1:idx],col=cols[2],pch=pchs[2])
  lines(mse.de.mestdagh.max[miR,1:idx],col=cols[3],pch=pchs[3])
  if (length(tcks.1>0)) {
    points(tcks.1,rep(minY,length(tcks.1)),pch=21,bg="black")
  }
  if (length(tcks.2>0)) {
    points(tcks.2,rep(minY,length(tcks.2)),pch=24,bg="grey")
    
  }
  legend("topleft",c("UCNR","LOD","MNV"),fill=cols) #,pch=pchs
}

plot.bias<-function(minY=NULL,maxY=NULL,idx=ncol(tobit.beta.est)) {
  par(mfrow=c(1,1))
  cols<-c("green","red","blue")#colors()[sample(100:150,5)]
  pchs<-c(21,22,23)#sample(19:25,3)  
  avg.mse.tobit<-apply(mse.de.tobit.varL[,1:idx],2,function(x) mean((x[x!=0])))
  avg.mse.mestdagh.max<-apply(mse.de.mestdagh.max[,1:idx],2, function(x) mean((x[x!=0])))
  avg.mse.mestdagh.varL<-apply(mse.de.mestdagh.varL[,1:idx],2, function(x) mean((x[x!=0])))
  if (is.null(maxY)) {
    maxY<-max(c(avg.mse.tobit,avg.mse.mestdagh.max,avg.mse.mestdagh.varL))
  }
  if (is.null(minY)) {
    minY<-min(c(avg.mse.tobit,avg.mse.mestdagh.max,avg.mse.mestdagh.varL))
  }
  plot(avg.mse.tobit,col=cols[1],pch=pchs[1],type="l",ylim=c(minY,maxY),main=expression(paste(italic(E(hat(delta) - delta))," per inclusion of censored observations")),ylab=expression(paste(italic(E(hat(delta) - delta)))),xlab="#Censored obs in dataset",lty=1,lwd=3)
  lines(avg.mse.mestdagh.varL,col=cols[2],pch=pchs[2],lty=2,lwd=3)
  lines(avg.mse.mestdagh.max,col=cols[3],pch=pchs[3],lty=3,lwd=3)
  abline(h=0)
  legend("topleft",c("UCNR","LOD","MNV"),col=cols,lty=1:3,lwd=3) #,pch=pchs
  grid(NA,4,lwd=2)
}


plot.bias.censstep<-function(s=ncol(tobit.beta.est),minY=NULL,maxY=NULL) {
  par(mfrow=c(1,1))
  cols<-c("green","red","blue")#colors()[sample(100:150,5)]
  pchs<-c(21,22,23)#sample(19:25,3)  
  avg.mse.tobit<-apply(mse.de.tobit.varL[,1:idx],2,function(x) mean(abs(x[x!=0])))
  avg.mse.mestdagh.max<-apply(mse.de.mestdagh.max[,1:idx],2, function(x) mean(abs(x[x!=0])))
  avg.mse.mestdagh.varL<-apply(mse.de.mestdagh.varL[,1:idx],2, function(x) mean(abs(x[x!=0])))
  if (is.null(maxY)) {
    maxY<-max(c(avg.mse.tobit,avg.mse.mestdagh.max,avg.mse.mestdagh.varL))
  }
  if (is.null(minY)) {
    minY<-min(c(avg.mse.tobit,avg.mse.mestdagh.max,avg.mse.mestdagh.varL))
  }
  
  plot(avg.mse.tobit,col=cols[1],pch=pchs[1],type="l",ylim=c(minY,maxY),main="Mean absolute bias per amount of censored obs",ylab="Mean absolute bias",xlab="#Censored obs in dataset")
  lines(avg.mse.mestdagh.varL,col=cols[2],pch=pchs[2])
  lines(avg.mse.mestdagh.max,col=cols[3],pch=pchs[3])
  legend("topleft",c("UCNR","LOD","MNV"),col=cols,lty=1:3,lwd=3) #,pch=pchs
}

plot.var<-function(minY=NULL,maxY=NULL,idx=ncol(tobit.beta.est)) {
  par(mfrow=c(1,1))
  cols<-c("green","red","blue")#colors()[sample(100:150,5)]
  pchs<-c(21,22,23)#sample(19:25,3)  
  avg.mse.tobit<-apply(mse.de.tobit.varL[,1:idx],2,function(x) mean(x[x!=0]^2))
  avg.mse.mestdagh.max<-apply(mse.de.mestdagh.max[,1:idx],2, function(x) mean(x[x!=0]^2))
  avg.mse.mestdagh.varL<-apply(mse.de.mestdagh.varL[,1:idx],2, function(x) mean(x[x!=0]^2))
  if (is.null(maxY)) {
    maxY<-max(c(avg.mse.tobit,avg.mse.mestdagh.max,avg.mse.mestdagh.varL))
  }
  if (is.null(minY)) {
    minY<-min(c(avg.mse.tobit,avg.mse.mestdagh.max,avg.mse.mestdagh.varL))
  }
  
  plot(avg.mse.tobit,col=cols[1],pch=pchs[1],type="l",ylim=c(minY,maxY),main=expression(paste(italic(E),(hat(delta) - delta)^2," per inclusion of censored observations")),ylab=expression(paste(italic(E),(hat(delta) - delta)^2)),xlab="#Censored obs in dataset",lty=1,lwd=3)
  lines(avg.mse.mestdagh.varL,col=cols[2],pch=pchs[2],lty=2,lwd=3)
  lines(avg.mse.mestdagh.max,col=cols[3],pch=pchs[3],lty=3,lwd=3)

  
  legend("topleft",c("UCNR","LOD","MNV"),col=cols,lty=1:3,lwd=3) #,pch=pchs
  grid(NA,4,lwd=2)
}

plot.var.est.mirs<-function(mirs=1:nrow(de.tobit.est.varL),minY=NULL,maxY=NULL,idx=ncol(tobit.beta.est)) {
  par(mfrow=c(1,1))
  cols<-c("green","red","blue")#colors()[sample(100:150,5)]
  pchs<-c(21,22,23)#sample(19:25,3)  
  var.mse.tobit<-apply(de.tobit.est.varL[mirs,1:idx],2,var)
  var.mse.mestdagh.max<-apply(de.mestdagh.est.max[mirs,1:idx],2, var)
  var.mse.mestdagh.varL<-apply(de.mestdagh.est.varL[mirs,1:idx],2, var)
  if (is.null(maxY)) {
    maxY<-max(c(var.mse.tobit,var.mse.mestdagh.max,var.mse.mestdagh.varL))
  }
  if (is.null(minY)) {
    minY<-min(c(var.mse.tobit,var.mse.mestdagh.max,var.mse.mestdagh.varL))
  }
    plot(var.mse.tobit,col=cols[1],pch=pchs[1],type="l",ylim=c(minY,maxY),main="Estimator variance per amount of censored obs",ylab="Variance",xlab="#Censored obs in dataset")
  lines(var.mse.mestdagh.varL,col=cols[2],pch=pchs[2])
  lines(var.mse.mestdagh.max,col=cols[3],pch=pchs[3])
  legend("topleft",c("UCNR","LOD","MNV"),fill=cols) #,pch=pchs
}

breakdown<-function(idx=ncol(tobit.beta.est),cens=T,diff=NULL) {
  which.cens<-apply(where.UV,1,function(x) any(x!=0))
  sum.cens<-apply(where.UV,1,function(x) sum(x==1|x==2))
  idx.cens<-which(which.cens==T)
  #    stop.pts<-apply(de.mestdagh.est.varL,1,function(x) which(x==0)[1])	
  bkd.tobit<-apply(conc.tobit[,1:idx],1, function(x) ifelse (any(diff(x)!=0),which(diff(x)!=0)[1],NA))
  bkd.mestdagh.varL<-apply(conc.mestdagh.varL[,1:idx],1, function(x) ifelse (any(diff(x)!=0),which(diff(x)!=0)[1],NA))
  bkd.mestdagh.max<-apply(conc.mestdagh.max[,1:idx],1, function(x) ifelse (any(diff(x)!=0),which(diff(x)!=0)[1],NA))
  
  stop.tobit<-rep(NA,length(bkd.tobit))
  stop.varL<-rep(NA,length(bkd.tobit))
  stop.max<-rep(NA,length(bkd.tobit))
  
  for (i in 1:length(bkd.tobit)) {
    if (!is.na(bkd.tobit[i])) {
      stop.tobit[i]<-(sum(where.UV[i,1:bkd.tobit[i]]>0))
    }
  }
  for (i in 1:length(bkd.mestdagh.varL)) {
    if (!is.na(bkd.mestdagh.varL[i])) {
      stop.varL[i]<-(sum(where.UV[i,1:bkd.mestdagh.varL[i]]>0))
    }
  }
  for (i in 1:length(bkd.mestdagh.max)) {
    if (!is.na(bkd.mestdagh.max[i])) {
      stop.max[i]<-(sum(where.UV[i,1:bkd.mestdagh.max[i]]>0))
    }
  }
  
  bkds<-cbind(stop.tobit/nrow(norm.est.ct),stop.varL/nrow(norm.est.ct),stop.max/nrow(norm.est.ct),sum.cens)
  colnames(bkds)<-c("UCNR","LOD","MNV","# Censored")
  if (cens) {
    print(bkds[idx.cens,])
  }
  if (!is.null(diff)) {
    print(bkds[diff,])
  }
}

plot.fdr<-function(alpha,nMi,idx=ncol(tobit.beta.est)) {
  fdr.alpha.max<-apply(p.mestdagh.est.max[,1:idx],2,function(x) sum(x<alpha))
  fdr.alpha.varL<-apply(p.mestdagh.est.varL[,1:idx],2,function(x) sum(x<alpha))
  fdr.alpha.tobit<-apply(p.tobit.est.varL[,1:idx],2,function(x) sum(x<alpha))
  plot(fdr.alpha.tobit/nMi,type="l",col="green",ylim=c(0,1))
  lines(fdr.alpha.varL/nMi,type="l",col="red")
  lines(fdr.alpha.max/nMi,type="l",col="blue")
}	

plot.pvalues<-function(miR,minY=-3,maxY=3,idx=ncol(de.tobit.est),minP=NULL,maxP=NULL,UCNR=T,UCNRW=F,MNV=F,MNV1=F,LOD=F,MOD=F,KNN=F,MODCT=F) {
  nf <- layout(matrix(c(1,1,2,3), 2, 2, byrow=TRUE), respect=TRUE)
  analyses<-c(UCNR,UCNRW,MNV,MNV1,LOD,MOD,KNN)
  cols<-c("green","red","blue","black","grey","yellow")
  buildLegend<-"UCNR"
  all.ps<-p.tobit.est.varL[miR,1:idx]
  des<-de.tobit.est.varL[miR,1:idx]
  #cols<-c("green","red","blue","black","grey","yellow")#colors()[sample(100:150,5)]"gray","wheat",
  plot(de.tobit.est.varL[miR,1:idx],col=cols[1],ylim=c(minY,maxY),main=paste("miRNA",miR," (final # UV : ",sum(where.UV[miR,]==1),"(MNA) /",sum(where.UV[miR,]==2)," (MNSC))"),ylab=expression(paste("Estimates Differential Expression ",hat(delta)[i])),xlab="#Censored obs in dataset",type="l",lwd=3,lty=1)
  howMany<-1
  if (UCNRW) {
    lines(de.tobit.est.weighted[miR,1:idx],col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"UCNRW")
    all.ps<-c(all.ps,p.tobit.est.weighted)
    des<-c(des,de.tobit.est.weighted[miR,1:idx])
  }
  if (MNV) {
    lines(de.mestdagh.est.max[miR,1:idx],col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"MNV")
    all.ps<-c(all.ps,p.mestdagh.est.max)
    des<-c(des,de.mestdagh.est.max[miR,1:idx])
  }
  if (MNV1) {
    lines(de.mestdagh.est.max1[miR,1:idx],col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"MNV+1")    
    all.ps<-c(all.ps,p.mestdagh.est.max1)
    des<-c(des,de.mestdagh.est.max1[miR,1:idx])
  }
  if (LOD) {
    lines(de.mestdagh.est.varL[miR,1:idx],col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"LOD")
    all.ps<-c(all.ps,p.mestdagh.est.varL)  
    des<-c(des,de.mestdagh.est.varL[miR,1:idx])
  }
  if (MOD) {
    lines(de.mestdagh.mod.glob[miR,1:idx],col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"MOD") 
    all.ps<-c(all.ps,p.mestdagh.mod.glob)
    des<-c(des,de.mestdagh.mod.glob[miR,1:idx])
  }
  if (MODCT) {
    lines(de.mestdagh.mod.glob.ct[miR,1:idx],col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"MODCT") 
    all.ps<-c(all.ps,p.mestdagh.mod.glob.ct)
    des<-c(des,de.mestdagh.mod.glob.ct[miR,1:idx])
  }
  if (KNN) {
    lines(de.mestdagh.est.knn[miR,1:idx],col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"KNN") 
    all.ps<-c(all.ps,p.mestdagh.est.knn)
    des<-c(des,de.mestdagh.est.knn[miR,1:idx])
  }
  abline(h=de.tobit.orig[miR])
  legend("topleft",buildLegend,col=cols[1:howMany],lty=1:howMany,lwd=3)
  grid(NA,4,lwd=2)
  tcks.1<-which(where.UV[miR,]==1)
  tcks.2<-which(where.UV[miR,]==2)
  if (length(tcks.1>0)) {
    points(tcks.1,rep(minY,length(tcks.1)),pch=21,bg="black")
  }
  if (length(tcks.2>0)) {
    points(tcks.2,rep(minY,length(tcks.2)),pch=24,bg="grey")
  }
  minY<-0
  maxY<-1
  if (is.null(minP)) {
    minP<-min(all.ps)-0.01
  }
  if (is.null(maxP)) {
    maxP<-1
  }
  plot(p.tobit.est.varL[miR,1:idx],type="l",col=cols[1],ylim=c(minP,maxP),main=paste("miRNA",miR," (# UV : ",sum(where.UV[miR,]==1),"/",sum(where.UV[miR,]==2),")"),ylab="P-values Differential Expression Hypothesis",xlab="#Censored obs in dataset",lty=1,lwd=3)
  howManyP<-1
  if (UCNRW) {
    lines(p.tobit.est.weighted[miR,1:idx],col=cols[howManyP+1],lwd=3,lty=howManyP+1)
    howManyP<-howManyP+1
  }
  if (MNV) {
    lines(p.mestdagh.est.max[miR,1:idx],col=cols[howManyP+1],lwd=3,lty=howManyP+1)
    howManyP<-howManyP+1
  }
  if (MNV1) {
    lines(p.mestdagh.est.max1[miR,1:idx],col=cols[howManyP+1],lwd=3,lty=howManyP+1)
    howManyP<-howManyP+1
  }
  if (LOD) {
    lines(p.mestdagh.est.varL[miR,1:idx],col=cols[howManyP+1],lwd=3,lty=howManyP+1)
    howManyP<-howManyP+1
  }
  if (MOD) {
    lines(p.mestdagh.mod.glob[miR,1:idx],col=cols[howManyP+1],lwd=3,lty=howManyP+1)
    howManyP<-howManyP+1
  }
  if (MODCT) {
    lines(p.mestdagh.mod.glob.ct[miR,1:idx],col=cols[howManyP+1],lwd=3,lty=howManyP+1)
    howManyP<-howManyP+1
  }
  if (KNN) {
    lines(p.mestdagh.est.knn[miR,1:idx],col=cols[howManyP+1],lwd=3,lty=howManyP+1)
    howManyP<-howManyP+1
  }
  
  nfs<-gl(howMany,idx,labels=buildLegend)
  boxplot(des~nfs,main=expression(paste("Box plot of ",hat(delta)[i])))
  abline(h=de.tobit.orig[miR])
}


plot.bias.all<-function(minY=NULL,maxY=NULL,idx=ncol(tobit.beta.est),UCNR=T,UCNRW=F,MNV=F,MNV1=F,LOD=F,MOD=F,KNN=F,MODCT=F) {
  par(mfrow=c(1,1))
  analyses<-c(UCNR,UCNRW,MNV,MNV1,LOD,MOD,MODCT,KNN)
  palette(rainbow(7))
  cols<-1:7
  buildLegend<-"UCNR"
  avg.mse.tobit<-apply(mse.de.tobit.varL[,1:idx],2,function(x) mean((x[x!=0])))
  avg.mse.mestdagh.max<-apply(mse.de.mestdagh.max[,1:idx],2, function(x) mean((x[x!=0])))
  avg.mse.mestdagh.max1<-apply(mse.de.mestdagh.max1[,1:idx],2, function(x) mean((x[x!=0])))
  avg.mse.mestdagh.varL<-apply(mse.de.mestdagh.varL[,1:idx],2, function(x) mean((x[x!=0])))
  avg.mse.mestdagh.knn<-apply(mse.de.mestdagh.knn[,1:idx],2, function(x) mean((x[x!=0])))
  avg.mse.tobit.weighted<-apply(mse.de.tobit.weighted[,1:idx],2, function(x) mean((x[x!=0])))
  avg.mse.mestdagh.mod.glob<-apply(mse.de.mestdagh.mod.glob[,1:idx],2, function(x) mean((x[x!=0])))
  avg.mse.mestdagh.mod.glob.ct<-apply(mse.de.mestdagh.mod.glob.ct[,1:idx],2, function(x) mean((x[x!=0])))
  if (is.null(maxY)) {
    maxY<-max(c(avg.mse.tobit,avg.mse.mestdagh.max,avg.mse.mestdagh.varL,avg.mse.mestdagh.knn,avg.mse.tobit.weighted,avg.mse.mestdagh.mod.glob,avg.mse.mestdagh.mod.glob.ct))
  }
  if (is.null(minY)) {
    minY<-min(c(avg.mse.tobit,avg.mse.mestdagh.max,avg.mse.mestdagh.varL,avg.mse.mestdagh.knn,avg.mse.tobit.weighted,avg.mse.mestdagh.mod.glob,avg.mse.mestdagh.mod.glob.ct))
  }
  plot(avg.mse.tobit,col="green",type="l",ylim=c(minY,maxY),main="",ylab=expression(paste(italic(E(hat(delta) - delta)))),xlab="#Censored obs in dataset",lty=1,lwd=3)
  howMany<-1
  if (UCNRW) {
    lines(avg.mse.tobit.weighted,col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"UCNRW")
  }
  if (MNV) {
    lines(avg.mse.mestdagh.max,col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"MNV")
  }
  if (MNV1) {
    lines(avg.mse.mestdagh.max1,col="blue",lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"MNV+1")    
  }
  if (LOD) {
    lines(avg.mse.mestdagh.varL,col="red",lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"LOD")
  }
  if (MOD) {
    lines(avg.mse.mestdagh.mod.glob,col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"MOD") 
  }
  if (MODCT) {
    lines(avg.mse.mestdagh.mod.glob.ct,col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"MODCT") 
  }
  if (KNN) {
    lines(avg.mse.mestdagh.knn,col="black",lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"KNN") 
  }
  abline(h=0)
  legend("topleft",buildLegend,col=c("green","blue","red","black"),lwd=3) #,pch=pchs
  #grid(NA,4,lwd=2)
}

plot.var.all<-function(minY=NULL,maxY=NULL,idx=ncol(tobit.beta.est),UCNR=T,UCNRW=F,MNV=F,MNV1=F,LOD=F,MOD=F,KNN=F,MODCT=F) {
  par(mfrow=c(1,1))
  analyses<-c(UCNR,UCNRW,MNV,MNV1,LOD,MOD,MODCT,KNN)
  palette(rainbow(7))
  cols<-1:7
  buildLegend<-"UCNR"
  avg.mse.tobit<-sqrt(apply(mse.de.tobit.varL[,1:idx],2,function(x) mean(x[x!=0]^2)))
  avg.mse.mestdagh.max<-sqrt(apply(mse.de.mestdagh.max[,1:idx],2, function(x) mean(x[x!=0]^2)))
  avg.mse.mestdagh.max1<-sqrt(apply(mse.de.mestdagh.max1[,1:idx],2, function(x) mean(x[x!=0]^2)))
  avg.mse.mestdagh.varL<-sqrt(apply(mse.de.mestdagh.varL[,1:idx],2, function(x) mean(x[x!=0]^2)))
  avg.mse.mestdagh.knn<-sqrt(apply(mse.de.mestdagh.knn[,1:idx],2, function(x) mean((x[x!=0]^2))))
  avg.mse.tobit.weighted<-sqrt(apply(mse.de.tobit.weighted[,1:idx],2, function(x) mean((x[x!=0]^2))))
  avg.mse.mestdagh.mod.glob<-sqrt(apply(mse.de.mestdagh.mod.glob[,1:idx],2, function(x) mean((x[x!=0]^2))))
  avg.mse.mestdagh.mod.glob.ct<-sqrt(apply(mse.de.mestdagh.mod.glob.ct[,1:idx],2, function(x) mean((x[x!=0]^2))))
  if (is.null(maxY)) {
    maxY<-max(c(avg.mse.tobit,avg.mse.mestdagh.max,avg.mse.mestdagh.varL,avg.mse.mestdagh.knn,avg.mse.tobit.weighted,avg.mse.mestdagh.mod.glob))
  }
  if (is.null(minY)) {
    minY<-min(c(avg.mse.tobit,avg.mse.mestdagh.max,avg.mse.mestdagh.varL,avg.mse.mestdagh.knn,avg.mse.tobit.weighted,avg.mse.mestdagh.mod.glob))
  }
    plot(avg.mse.tobit,col="green",type="l",ylim=c(minY,maxY),main="",ylab=expression(sqrt(paste(italic(E),(hat(delta) - delta)^2))),xlab="#Censored obs in dataset",lty=1,lwd=3)
  howMany<-1
  if (UCNRW) {
    lines(avg.mse.tobit.weighted,col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"UCNRW")
  }
  if (MNV) {
    lines(avg.mse.mestdagh.max,col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"MNV")
  }
  if (MNV1) {
    lines(avg.mse.mestdagh.max1,col="blue",lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"MNV+1")    
  }
  if (LOD) {
    lines(avg.mse.mestdagh.varL,col="red",lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"LOD")
  }
  if (MOD) {
    lines(avg.mse.mestdagh.mod.glob,col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"MOD") 
  }
  if (MODCT) {
    lines(avg.mse.mestdagh.mod.glob.ct,col=cols[howMany+1],lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"MODCT") 
  }
  if (KNN) {
    lines(avg.mse.mestdagh.knn,col="black",lwd=3,lty=howMany+1)
    howMany<-howMany+1
    buildLegend<-c(  buildLegend,"KNN") 
  }
  
  
  legend("topleft",buildLegend,col=c("green","blue","red","black"),lty=1:howMany,lwd=3) #,pch=pchs
  #grid(NA,4,lwd=2)
}