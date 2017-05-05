### model matrix for UCNR
#########################

unified.design.matrix<-function(miRNA,subject,group,nrMi,nrS,intercept=TRUE,sub.sum=F,UV=NULL) {
  contrasts(miRNA)<-contr.sum
  # sum coding subject
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



### censored normal heteroscedastic regression (function)
#########################################################

censHeteroLogLikSigma <- function(param) {
  mu <- param[1:ncol(X)]
  sigma <- param[(ncol(X)+1):(ncol(X)+ncol(Z))]
  yHat<-X%*%mu
  sigmaHat<-exp((Z%*%sigma))
  ll <- rep(NA, length(cq))
  ll[obsBelow] <- dnorm((cq - yHat)[obsBelow]/sigmaHat[obsBelow], log = TRUE) - log(sigmaHat)[obsBelow]
  ll[obsAbove] <- pnorm((yHat[obsAbove] - right)/sigmaHat[obsAbove], log.p = TRUE)
  
  grad <- matrix(NA, nrow = length(cq), ncol = length(param))
  grad[obsBelow, ] <- cbind(((cq - yHat)[obsBelow]/sigmaHat[obsBelow]) * X[obsBelow, , drop = FALSE]/sigmaHat[obsBelow], (((cq - yHat)[obsBelow]/sigmaHat[obsBelow])^2 - 1)* Z[obsBelow, , drop = FALSE] /2)
  grad[obsAbove, ] <- dnorm((yHat[obsAbove] - right)/sigmaHat[obsAbove])/pnorm((yHat[obsAbove] - right)/sigmaHat[obsAbove]) * cbind(X[obsAbove, , drop = FALSE]/sigmaHat[obsAbove], (-(yHat[obsAbove] - right)/sigmaHat[obsAbove])* cbind(Z[obsAbove, , drop = FALSE]))
  
  attr(ll, "gradient") <- grad
  return(ll)
}

### UCNR 
#############################
ucnr.fit <- function(data, truncated){
	### set LOD at truncated
	cat("Preparing for model fitting...\n")
	cq<-data$Cq;
	cq<<-replace(cq,cq>=truncated,truncated)
	
	### extract necessary vectors
	group<-factor(data$group)
	genes<-factor(data$target)
	subject<-factor(data$sample)
	mm<-unified.design.matrix(genes,subject,group,nlevels(genes),nlevels(subject),sub.sum=F)
	
	cat("Fitting UCNR model...\n")
	#### Fit Tobit model for starting values
	tobit.fitA<-tobit(cq~mm-1,left=-Inf,right=truncated)

	### Fit UCNR
	right<<-truncated
	obsAbove <<- cq >= right
	obsBelow <<- !obsAbove
	X<<-mm
	Z<<-model.matrix(~genes,contrasts=list(genes="contr.treatment"))
	param<-c(coef(tobit.fitA),c(tobit.fitA$scale,rep(0,ncol(Z)-1)))
	result <- maxLik(censHeteroLogLikSigma, start=c(coef(tobit.fitA),c(tobit.fitA$scale,rep(0,ncol(Z)-1))),method="BFGS")
	
	cat("Processing UCNR model...\n")
	estimates <- summary(result)$estimate[(nlevels(genes)+nlevels(subject)):(nlevels(genes)+nlevels(subject)+nlevels(genes)-2),]
	return(estimates)
}	

### UCNR with reference genes
#############################

ucnr.fit.refgenes <- function(data, mm, truncated){
	cat("Preparing for model fitting...\n")
	genes<- data$target
	subject<-data$sample
	array<-data$sample
	refG<-data$referenceGeneID
	group<-as.factor(data$group)
	
	group1<- subject[group=="0"]
	group2<- subject[group=="1"]
	
	matR<- matrix(data$referenceGeneID,nc=nlevels(genes))
	refId<-which(apply(matR,2,function (x) sum(x))!=0)
	nrOfNonRefGenes<-nlevels(genes)-length(refId)
	nonRefId<-1:nlevels(genes)
	nonRefId<-nonRefId[-refId]
	

	# remove last gene (due to contrast coding)
	toDel<- c(1 + nlevels(genes)-1 + nlevels(subject)-1 + nlevels(array)-1 + nrOfNonRefGenes)

	toImpStart<-toDel-nrOfNonRefGenes+1
	toImpStop<-toDel-1

	mm[,toImpStart:toImpStop]<- mm[,toImpStart:toImpStop]-mm[,toDel]
	mm<-mm[,-toDel]


	genes.resp<-c()  
	for (i in 1:length(nonRefId)) {
		genes.resp<-c(genes.resp,((nonRefId[i]-1)*nlevels(subject)+1):(nonRefId[i]*nlevels(subject)))
	}
	mm[genes.resp,(1 + nlevels(genes)-1 + nlevels(subject)):(1 + nlevels(genes)-1 + nlevels(subject)-1 + nlevels(array)-1)]<-0


	data$Cq<-replace(data$Cq,is.na(data$Cq),truncated)
	data$Cq<-replace(data$Cq,data$Cq>truncated,truncated)
	data$cq<-data$Cq
	cat("Fitting UCNR model...\n")
	### fit tobit model for starting values
	X<<-mm
	Z<<-model.matrix(~target, data=data) 
	fit.genes<-tobit(cq~X-1,data=data,left=-Inf,right=truncated)
	coef.tobit<-coef(summary(fit.genes))[(nlevels(genes)+nlevels(subject)+nlevels(array)-1):(nlevels(genes)+nlevels(subject)+nlevels(array)-1+nlevels(genes)-7),]

	# maximum likelihood estimation for censored normal heteroscedastic regression
	require(maxLik)
	# set some parms
	cq<<-data$Cq
	right<<-truncated
	obsAbove <<- cq >= right
	obsBelow <<- !obsAbove
	# maxLik
	
	maxLikCHLLS <- maxLik(censHeteroLogLikSigma, start=c(coef(fit.genes),c(log(fit.genes$scale),rep(0,ncol(Z)-1))),method="BFGS")



	cat("Processing UCNR model...\n")
	cfs.table.chr<-suppressWarnings(coef(summary(maxLikCHLLS)))
	mu<-cfs.table.chr[1,1]
	alpha<-cfs.table.chr[2:nlevels(genes),1]
	alpha<-c(alpha,-sum(alpha))
	
	idx.eta.1<-c()
	for (a in 1:15) {
		idx.eta.1[a]<-grep(unique(group1)[a],rownames(cfs.table.chr))[2]
	}
	eta1<-cfs.table.chr[idx.eta.1,1]
	eta1<-replace(eta1,is.na(eta1),0)
	
	idx.eta.2<-c()
	for (a in 1:15) {
	  idx.eta.2[a]<-grep(unique(group2)[a],rownames(cfs.table.chr))[2]
	}
	eta2<-cfs.table.chr[idx.eta.2,1]
	eta2<-replace(eta2,is.na(eta2),0)

	delta<-cfs.table.chr[(nlevels(genes)+2*nlevels(subject)-1):(nrow(cfs.table.chr)-nlevels(genes)),1]
	delta<-c(delta,-sum(delta))
	delta_i<-delta+mean(eta1)-mean(eta2)     

	## construct contrast matrix
	B<-cfs.table.chr[1:(nrow(cfs.table.chr)-nlevels(genes)),1]
	covB<--solve(maxLikCHLLS$hessian)[1:(nrow(cfs.table.chr)-1),1:(nrow(cfs.table.chr)-1)]
	covB<-covB[1:ncol(X),1:ncol(X)] # fix dim
	dim(covB)
	length(B)
	L<-rep(0,length(B))#matrix(0,nc=length(B),nr=31)
	
	B[idx.eta.1]
	B[idx.eta.2]
	L[idx.eta.1[2:15]]<-1/15
	L[idx.eta.2]<--1/15
	
	# SE
	sqrt(L%*%covB%*%L)
	
	# Wald Test
	(L%*%B)*solve(L%*%covB%*%L)*(L%*%B)
	
	delta_L<-c()
	SE<-c()
	WS<-c()
	pvalues<-c()
	idx.delta<-which(names(B)=="XgenesAHCY:group1")

	for (i in 1:(length(delta))-1) {
		L<-rep(0,length(B))
		L[idx.eta.1[2:15]]<-1/15
		L[idx.eta.2]<--1/15
		L[idx.delta+i-1]<-1
		delta_L[i]<-L%*%B
		# SE
		SE[i]<-sqrt(L%*%covB%*%L)
		# Wald Test
		WS[i]<-(L%*%B)*solve(L%*%covB%*%L)*(L%*%B)
		pvalues[i]<-1-pchisq(WS[i],1)
	}

	coef.contrasts<-cbind(delta_L,SE,WS,pvalues)
	colnames(coef.contrasts)<-c("Estimate","SE","Wald Stat","p value")
	rownames(coef.contrasts)<-names(B)[idx.delta:178]

	L<-rep(0,length(B))
	L[idx.eta.1[2:15]]<-1/15
	L[idx.eta.2]<--1/15
	L[(idx.delta):(idx.delta+(length(delta)-1)-1)]<--1
	delta_L[length(delta)]<-L%*%B
	# SE
	SE[length(delta)]<-sqrt(L%*%covB%*%L)
	# Wald Test
	WS[length(delta)]<-(L%*%B)*solve(L%*%covB%*%L)*(L%*%B)
	pvalues[length(delta)]<-1-pchisq(WS[length(delta)],1)
	
	coef.contrasts<-cbind(delta_L,SE,WS,pvalues)
	colnames(coef.contrasts)<-c("Estimate","SE","Wald Stat","p value")
	rownames(coef.contrasts)<-names(B)[idx.delta:179]

	adj.pvalues<-mt.rawp2adjp(coef.contrasts[,4], proc=c("BH"))  
	adj.pvalues<-adj.pvalues$adjp[order(adj.pvalues$index),2]
	coef.contrasts<-cbind(coef.contrasts,adj.pvalues)
	rownames(coef.contrasts)<-as.character(unique(genes)[nonRefId])	
	return(as.data.frame(coef.contrasts))
}