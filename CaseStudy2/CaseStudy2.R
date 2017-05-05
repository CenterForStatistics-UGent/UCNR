### Case study 2
#######################################

require(AER)
require(maxLik)
require(multtest)

Risk <- read.csv("data/SIOPEN.csv")
source("scripts/ucnrfunctions.R")

### order dataset; first target genes and then reference genes, then alphabetical order
Risk<-Risk[order(Risk$referenceGeneID,Risk$target),]

### rename, the dataset needs following columns:
###### sample: sample identifier
###### target: gene identifier
###### Cq: quantification cycle
###### group: two groups to compare, coded as 0/1
###### referenceGeneID: whether the gene is a reference gene or a target gene
colnames(Risk)[13] <- "group"

### Differential Expression Analysis using UCNR
###############################################
genes<-Risk$target
subject<-Risk$sample
array<-Risk$sample
refG<-Risk$referenceGeneID
group<-as.factor(Risk$group)
contrasts(genes)<-contr.sum

### Construct unified model matrix
mm<-model.matrix(~genes+subject+array+genes:group)

matR<- matrix(Risk$referenceGeneID,nc=nlevels(genes))
refId<-which(apply(matR,2,function (x) sum(x))!=0)
nrOfNonRefGenes<-nlevels(genes)-length(refId)
nonRefId<-1:nlevels(genes)
nonRefId<-nonRefId[-refId]

### Which are the reference genes? 
names(sort(table(Risk[Risk$referenceGeneID==1,3]),decreasing=T)[1:5])

### remove interaction effects for reference genes
genes.toRemove<- which(colnames(mm)%in%c("genesAluSq:group1","genesHMBS:group1","genesHPRT1:group1","genesSDHA:group1","genesUBC:group1"))
mm<-mm[,-genes.toRemove]

### set LOD at 39
truncated<-39

### fit the UCNR model with reference genes
fit <- ucnr.fit.refgenes(data=Risk, mm=mm, truncated=truncated)
fit

### Differential Expression Analysis using LOD/MNV+1
####################################################
Risk.mat<-matrix(Risk$Cq,nc=63)
colnames(Risk.mat)<-unique(Risk$target)
rownames(Risk.mat)<-Risk$sample[1:30]

### Mean Expression Value
MEN<-apply(Risk.mat[,59:63],1,function(x) mean(x[x<truncated]))
### Normalize according to LOD
norm.LOD<-apply(Risk.mat[,1:58],2,function(x) x-MEN)
### Normalize according to MNV
norm.MNV<-replace(Risk.mat,Risk.mat>=truncated,NA)
norm.MNV<-apply(norm.MNV[,1:58],2,function(x) x-MEN)
max.mi<-apply(norm.MNV,2,function(x) max(x,na.rm=T))

for (i in 1:ncol(norm.MNV)) {
  
  ident.na<-which(is.na(norm.MNV[,i]))
  if (length(ident.na) > 0) {
    norm.MNV[,i]<-replace(norm.MNV[,i],is.na(norm.MNV[,i]),max.mi[i])
  }    
}

### Prepare for t-test
checkInf<-apply(norm.MNV,2,function(x) any(x==-Inf))
idNotToUse<-which(checkInf)
checkEq<-apply(norm.MNV,2,function(x) sum(x[1:15]!=x[16:30]))
idAllEqual<-which(checkEq==0)
### Perform t-test
ttestlist.LOD<-apply(norm.LOD,2,function(x) {t.test(x[1:15],x[16:30])})
ttestlist.MNV<-apply(norm.MNV,2,function(x) {t.test(x[1:15],x[16:30])})

length(ttestlist.LOD)
length(ttestlist.MNV)

### Extract p-values and DE
pvalues.LOD<-lapply(ttestlist.LOD,function(x) x$p.value)
pvalues.MNV<-lapply(ttestlist.MNV,function(x) x$p.value)
de.LOD<-lapply(ttestlist.LOD,function(x) diff(x$estimate))
de.MNV<-lapply(ttestlist.MNV,function(x) diff(x$estimate))

### combine
LODoutput<-cbind(unlist(de.LOD),unlist(pvalues.LOD))
MNVoutput<-cbind(unlist(de.MNV),unlist(pvalues.MNV))

### calculate adjuested p-values
adj.pvalues.LOD<-mt.rawp2adjp(LODoutput[,2], proc=c("BH"))  
adj.pvalues.LOD<-adj.pvalues.LOD$adjp[order(adj.pvalues.LOD$index),2]
LODoutput<-cbind(LODoutput,adj.pvalues.LOD)

adj.pvalues.MNV<-mt.rawp2adjp(MNVoutput[,2], proc=c("BH"))  
adj.pvalues.MNV<-adj.pvalues.MNV$adjp[order(adj.pvalues.MNV$index),2]
MNVoutput<-cbind(MNVoutput,adj.pvalues.MNV)

rownames(LODoutput)<-names(unlist(pvalues.LOD))
rownames(MNVoutput)<-names(unlist(pvalues.MNV))

### number of significant genes (adjusted p-values < 0.05)
dim(fit[fit[,5]<0.05,]) 
dim(LODoutput[LODoutput[,3]<0.05,])
dim(MNVoutput[MNVoutput[,3]<0.05,])