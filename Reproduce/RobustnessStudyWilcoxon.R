###############################################################
### Script for consistency study - miRs 
###			Each iteration 1 miR is censored
###############################################################

#make sure to set correct locale for exact reproduction
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")

require(samr)
require(multtest)
require(survival)
require(e1071)
require(lmtest)
require(coin)
require("lattice")
require(brglm)
require(AER)
require(grid)
require(SLqPCR)
require(impute)

path<-"";
srcpath<-paste(path,"data/",sep="");
outpath<-paste(path,"results/",sep="");
scriptpath<-paste(path,"scripts/",sep="");

source(paste(scriptpath,"read in neuroblastoma.R",sep=""))
source(paste(scriptpath,"auxfunctions.R",sep=""))

### the number of miRNAs in the dataset: set to 50 for main paper results
### change to 100/200 for appendix results to demonstrate scalability
testNrMis<-50


#############
# Read data #
#############

names(NB.ok.stack)[1]<-"subject";
ok.stack<-NB.ok.stack

# Consider group MYCNSC
NB.mat.mycnsc<-NB.mat[mycnsc,]
# Mean Expression Value
NB.mean.expr<-apply(NB.mat.mycnsc,1,function(x) mean(x[x!=40]))
# Normalize
NB.norm<-apply(NB.mat.mycnsc,2,function(x) x-NB.mean.expr)
# Count #UD each miRNA
UD.miRNA<-apply(NB.mat.mycnsc,2,function(x) sum(ifelse(x==40,1,0)))

#which miRNAs have no UV? call this all expressed miRNAs
all.expr.miRNAs<-names(which(UD.miRNA==0))
#get row IDs for all-expressed miRNAs
idx.all<-c()
for (i in 1:length(all.expr.miRNAs)) {
  idx.all<-c(idx.all,which(NB.stack$miRNA==all.expr.miRNAs[i]))
}
#get all features for these all-expr miRNAs
NB.stack.all.expr<-NB.stack[idx.all,]

idx.over<-which(NB.stack.all.expr[,5]==1)
problem.miRNAs<-unique(NB.stack.all.expr[idx.over,2])
idx.not<-c()
for (i in 1:length(problem.miRNAs)) {
  idx.not<-c(idx.not,which(NB.stack.all.expr$miRNA==problem.miRNAs[i]))
}
NB.stack.all.expr<-NB.stack.all.expr[-idx.not,]
NB.stack.all.expr<-NB.stack.all.expr[order(NB.stack.all.expr$miRNA,NB.stack.all.expr$group),]




#extract the set number of miRNAs,
testData<-NB.stack.all.expr[1:(61*testNrMis),]

# Delta and nr Diff
delta.cq<-2
nDiff<-10
S<-1100
alpha<-0.05

group<-factor(testData$group)
miRNA<-factor(rep(1:length(unique(testData$miRNA)),each=length(unique(testData$NB))))
subject<-factor(rep(1:length(unique(testData$NB)),length(unique(testData$miRNA))))

nMi<-length(unique(miRNA))
nS<-length(unique(subject))

cq.grp1<-matrix(testData[,3],nc=nMi)[1:22,]
cq.grp2<-matrix(testData[,3],nc=nMi)[23:61,]

# Compute shift
mean.cq.grp1<-apply(cq.grp1,2,mean)
mean.cq.grp2<-apply(cq.grp2,2,mean)
shift<-mean.cq.grp1-mean.cq.grp2

for(i in 1:nMi) {  
  cq.grp2[,i]<-cq.grp2[,i]+shift[i]
}
# store in testData$Ct
joined.grps<-rbind(cq.grp1,cq.grp2)
testData$Ct<-c(joined.grps)


# Normalisation Constant
norm.ct<-rowMeans(matrix(testData$Ct,nc=nMi))

# Normalise data
normData<-matrix(testData$Ct,nc=nMi)-norm.ct
norm.cq.grp1<-normData[1:22,]
norm.cq.grp2<-normData[23:61,]


# Add Delta to normalised data
for(i in 1:nMi) {	
  #	norm.cq.grp2[,i]<-norm.cq.grp2[,i]+shift[i]
  if (i <=(nDiff)) { 
    norm.cq.grp2[,i]<-norm.cq.grp2[,i]+delta.cq
  }
  else if (i <=2*nDiff) {
    norm.cq.grp2[,i]<-norm.cq.grp2[,i]-delta.cq
  }
}

# transform back to use Cq

joined.grps<-rbind(norm.cq.grp1,norm.cq.grp2)
joined.grps<-joined.grps+norm.ct
testData$Ct<-c(joined.grps)

# or analysis :
LOD<-40

contrasts(testData$NB)<-contr.sum
tobit.subj<-coef(tobit(Ct~NB,data=testData))
tobit.norm.ct<-c()
for (i in 2:length((unique(testData$NB)))) {
  tobit.norm.ct[i-1]<-tobit.subj[1]+tobit.subj[i]
}
# similar as
norm.ct<-rowMeans(matrix(testData$Ct,nc=nMi))

mm.unified<-unified.design.matrix(miRNA,subject,group,nlevels(miRNA),nlevels(subject),sub.sum=T)
tobit.subj.unified<-coef(tobit(Ct~mm.unified-1,data=testData,left=-Inf,right=LOD))

tobit.beta<-c()
for (i in 2:nS) {
  tobit.beta[i-1]<-tobit.subj.unified[1]+tobit.subj.unified[i+nMi-1]
}

mm.mestdagh<-mestdagh.design.matrix(miRNA,group,nrMi=nMi,nrS=nS,intercept=TRUE,miRNA.sum=T)



norm.est.ct<-apply(matrix(testData$Ct,nc=nMi),1,function(x) mean(x[x<LOD]))
normData<-matrix(testData$Ct,nc=nMi)-norm.est.ct
de.mestdagh.orig<-c()
de.mestdagh.orig<-coef(lm(c(normData)~mm.mestdagh-1))[(nMi+1):(2*nMi-1)]
de.tobit.orig<-c()
de.tobit.orig<-tobit.subj.unified[(nMi+nS):(2*nMi+nS-2)]

orig.cq<-testData$Ct
cq.during.sim<-orig.cq
cq.during.sim.varL<-orig.cq

### 
# Impute s LODs

de.mestdagh.est<-matrix(0,nr=nMi-1,nc=S)
de.tobit.est<-matrix(0,nr=nMi-1,nc=S)
tobit.beta.est<-matrix(0,nr=nS-1,nc=S)
norm.est.ct<-matrix(0,nr=nS,nc=S)
norm.mod.glob.est<-matrix(0,nr=nS,nc=S)
norm.mod.glob.est.ct<-matrix(0,nr=nS,nc=S)
where.UV<-matrix(0,nr=nMi,nc=S)
mse.de.tobit<-c()
mse.de.mestdagh<-c()
mse.nc.tobit<-c()
mse.nc.mestdagh<-c()
de.mestdagh.est.varL<-matrix(0,nr=nMi,nc=S)
de.mestdagh.est.max<-matrix(0,nr=nMi,nc=S)
de.mestdagh.est.max1<-matrix(0,nr=nMi,nc=S)
de.mestdagh.est.knn<-matrix(0,nr=nMi,nc=S)
de.mestdagh.mod.glob<-matrix(0,nr=nMi,nc=S)
de.mestdagh.mod.glob.ct<-matrix(0,nr=nMi,nc=S)
p.mestdagh.est.knn<-matrix(1,nr=nMi,nc=S)
p.mestdagh.est.max<-matrix(1,nr=nMi,nc=S)
p.mestdagh.est.max1<-matrix(1,nr=nMi,nc=S)
p.mestdagh.est.varL<-matrix(1,nr=nMi,nc=S)
p.mestdagh.mod.glob<-matrix(1,nr=nMi,nc=S)
p.mestdagh.mod.glob.ct<-matrix(1,nr=nMi,nc=S)
ciL.mestdagh.est.varL<-matrix(0,nr=nMi,nc=S)
ciU.mestdagh.est.varL<-matrix(0,nr=nMi,nc=S)
ciL.mestdagh.est.max<-matrix(0,nr=nMi,nc=S)
ciU.mestdagh.est.max<-matrix(0,nr=nMi,nc=S)
std.mestdagh.est.varL<-matrix(0,nr=nMi,nc=S)
std.mestdagh.est.max<-matrix(0,nr=nMi,nc=S)
de.tobit.est.varL<-matrix(0,nr=nMi-1,nc=S)
de.tobit.est.weighted<-matrix(0,nr=nMi-1,nc=S)
p.tobit.est.varL<-matrix(1,nr=nMi-1,nc=S)
p.tobit.est.weighted<-matrix(1,nr=nMi-1,nc=S)
std.tobit.est.varL<-matrix(1,nr=nMi-1,nc=S)
tobit.beta.est.varL<-matrix(0,nr=nS-1,nc=S)
mse.de.tobit.varL<-matrix(0,nr=nMi-1,nc=S)
mse.de.tobit.weighted<-matrix(0,nr=nMi-1,nc=S)
mse.de.mestdagh.varL<-matrix(0,nr=nMi,nc=S)
mse.de.mestdagh.max<-matrix(0,nr=nMi,nc=S)
mse.de.mestdagh.max1<-matrix(0,nr=nMi,nc=S)
mse.de.mestdagh.knn<-matrix(0,nr=nMi,nc=S)
mse.de.mestdagh.mod.glob<-matrix(0,nr=nMi,nc=S)
mse.de.mestdagh.mod.glob.ct<-matrix(0,nr=nMi,nc=S)
mse.nc.tobit.varL<-matrix(0,nr=nS-1,nc=S)
mse.nc.mestdagh.varL<-matrix(0,nr=nS,nc=S)
conc.tobit<-matrix(0,nr=nMi,nc=S)
conc.tobit.weighted<-matrix(0,nr=nMi,nc=S)
conc.mestdagh.varL<-matrix(0,nr=nMi,nc=S)
conc.mestdagh.max<-matrix(0,nr=nMi,nc=S)
conc.mestdagh.max1<-matrix(0,nr=nMi,nc=S)
conc.mestdagh.knn<-matrix(0,nr=nMi,nc=S)
conc.mestdagh.mod.glob<-matrix(0,nr=nMi,nc=S)
conc.mestdagh.mod.glob.ct<-matrix(0,nr=nMi,nc=S)
ind.DE<-1:(nMi-1)
ind.DEt<-1:(nMi)
vLOD<-LOD
sink(file="Rlog_Wilcox50miRNAs.txt", split=TRUE)
for (s in 1:S) {
  
  miRNA.cens<-ceiling(which(cq.during.sim.varL==max(cq.during.sim.varL[cq.during.sim.varL!=vLOD]))/61)
  miRNA.grp.cens<-ifelse(which(cq.during.sim.varL==max(cq.during.sim.varL[cq.during.sim.varL!=vLOD]))%%61>22,2,1)
  where.UV[ind.DE[miRNA.cens],s]<-miRNA.grp.cens
  
  vLOD<-max(cq.during.sim.varL[cq.during.sim.varL!=vLOD]) 					# hlp is max, but less than LOD
  cq.during.sim.varL<-replace(cq.during.sim.varL,cq.during.sim.varL>=vLOD,NA)  # replace max < LOD
  
  print(paste(s,": ",LOD," - ",vLOD))
  UV.grp1<-apply(matrix(cq.during.sim.varL,nc=nMi)[1:22,],2,function(x) sum(is.na(x)))
  UV.grp2<-apply(matrix(cq.during.sim.varL,nc=nMi)[23:61,],2,function(x) sum(is.na(x)))
  UV.grp<-apply(matrix(cq.during.sim.varL,nc=nMi),2,function(x) sum(is.na(x)))
  print(paste(UV.grp))
  if(any(UV.grp1>=21) | any(UV.grp2>=38)) {
    idx.miRNA.1<-which(UV.grp1>=21)
    idx.miRNA.2<-which(UV.grp2>=38)
    if (length(idx.miRNA.1) > 0) {
      indices<-c()
      print(paste("MiRNA",idx.miRNA.1," in grp 1 is out",UV.grp1[idx.miRNA.1]))
      cq.during.sim.varL<-c(matrix(cq.during.sim.varL,nc=nMi)[,-idx.miRNA.1])
      for (ind in 1:length(idx.miRNA.1)){
        indices<-c(indices,((idx.miRNA.1-1)*61+1):(idx.miRNA.1*61))		
      }
      miRNA<-factor(miRNA[-indices])
      subject<-factor(subject[-indices])
      group<-factor(group[-indices])
      ind.DE<-ind.DE[-idx.miRNA.1]
      ind.DEt<-ind.DEt[-idx.miRNA.1]		
      print(paste("Indices left over",ind.DE))
    }
    else if (length(idx.miRNA.2) > 0) {
      indices<-c()
      print(paste("MiRNA",idx.miRNA.2," is grp 2 is out",UV.grp2[idx.miRNA.2]))
      cq.during.sim.varL<-c(matrix(cq.during.sim.varL,nc=nMi)[,-idx.miRNA.2])
      for (ind in 1:length(idx.miRNA.2)){
        indices<-c(indices,((idx.miRNA.2-1)*61+1):(idx.miRNA.2*61))
      }
      miRNA<-factor(miRNA[-indices])
      subject<-factor(subject[-indices])
      group<-factor(group[-indices])
      ind.DE<-ind.DE[-idx.miRNA.2]
      ind.DEt<-ind.DEt[-idx.miRNA.2]	
      print(paste("Indices left over",ind.DE))
    }
    nMi<-length(unique(miRNA))#
    mm.unified<-unified.design.matrix(miRNA,subject,group,nlevels(miRNA),nlevels(subject),sub.sum=T)
    mm.mestdagh<-mestdagh.design.matrix(miRNA,group,nrMi=nMi,nrS=nS,intercept=TRUE,miRNA.sum=T)
  }
  
  print("Voor knn\n")
  print("modified global mean\n")
  norm.mod.glob<-scale(matrix(cq.during.sim.varL,nc=nMi),center=T,scale=F)
  scaleFactors<-colMeans(matrix(cq.during.sim.varL,nc=nMi),na.rm=T)
  norm.mod.glob.est[,s]<-rowMeans(norm.mod.glob,na.rm=T)
  normData.mod.glob<-norm.mod.glob-norm.mod.glob.est[,s]
  cnt.na.miR<-apply(norm.mod.glob,2,function(x) sum(is.na(x)))
  only.expr<-which(cnt.na.miR==0)
  norm.mod.glob.est.ct[,s]<-rowMeans(norm.mod.glob[,only.expr])#rowMeans(matrix(testData$Ct,nc=nMi)) #
  normData.mod.glob.ct<-norm.mod.glob-norm.mod.glob.est.ct[,s]
  
  normData.varL<-normData.mod.glob
  normData.max<-normData.mod.glob
  normData.max1<-normData.mod.glob
  
  normData.knn<-t(normData.max)
  normData.knn<-impute.knn(normData.knn,k=10)$data
  normData.knn<-t(normData.knn)
  
  print("Na knn\n")
  
  # check if all UV in groups appear
  for (i in 1:nMi) {
    ident.na<-which(is.na(normData.varL[,i]))
    if (length(ident.na) > 0) {
    #	print(i)
      max.mi<-max(normData.varL[,i],na.rm=T)
      #print(max.mi)
      max.mi.mod.glob<-max(normData.mod.glob[,i],na.rm=T)
      #print(max.mi.mod.glob)
      max.mi.mod.glob.ct<-max(normData.mod.glob.ct[,i],na.rm=T)
      normData.varL[,i]<-replace(normData.varL[,i],is.na(normData.varL[,i]),vLOD-scaleFactors[i])
      normData.max[,i]<-replace(normData.max[,i],is.na(normData.max[,i]),max.mi)
      normData.max1[,i]<-replace(normData.max1[,i],is.na(normData.max1[,i]),max.mi+1)
      normData.mod.glob[,i]<-replace(normData.mod.glob[,i],is.na(normData.mod.glob[,i]),max.mi.mod.glob)
      #print(normData.mod.glob[,i])
      normData.mod.glob.ct[,i]<-replace(normData.mod.glob.ct[,i],is.na(normData.mod.glob.ct[,i]),max.mi.mod.glob.ct)      
    }		
  }
  
  print("Voor ttest\n")
  de.mestdagh.est.varL[ind.DEt,s]<-apply(normData.varL,2,function(x) diff(t.test(x[1:22],x[23:61])$estimate))
  #ciL.mestdagh.est.varL[ind.DEt,s]<-apply(normData.varL,2,function(x) t.test(x[23:61],x[1:22])$conf.int[1])
  #ciU.mestdagh.est.varL[ind.DEt,s]<-apply(normData.varL,2,function(x) t.test(x[23:61],x[1:22])$conf.int[2])
  p.mestdagh.est.varL[ind.DEt,s]<-apply(normData.varL,2,function(x) wilcox.test(x[1:22],x[23:61])$p.value)
  
  de.mestdagh.est.max[ind.DEt,s]<-apply(normData.max,2,function(x) diff(t.test(x[1:22],x[23:61])$estimate))
  #ciL.mestdagh.est.max[ind.DEt,s]<-apply(normData.max,2,function(x) t.test(x[23:61],x[1:22])$conf.int[1])
  #ciU.mestdagh.est.max[ind.DEt,s]<-apply(normData.max,2,function(x) t.test(x[23:61],x[1:22])$conf.int[2])
  p.mestdagh.est.max[ind.DEt,s]<-apply(normData.max,2,function(x) wilcox.test(x[1:22],x[23:61])$p.value)
  
  de.mestdagh.est.max1[ind.DEt,s]<-apply(normData.max1,2,function(x) diff(t.test(x[1:22],x[23:61])$estimate))
  #ciL.mestdagh.est.max[ind.DEt,s]<-apply(normData.max,2,function(x) t.test(x[23:61],x[1:22])$conf.int[1])
  #ciU.mestdagh.est.max[ind.DEt,s]<-apply(normData.max,2,function(x) t.test(x[23:61],x[1:22])$conf.int[2])
  p.mestdagh.est.max1[ind.DEt,s]<-apply(normData.max1,2,function(x) wilcox.test(x[1:22],x[23:61])$p.value)

  de.mestdagh.est.knn[ind.DEt,s]<-apply(normData.knn,2,function(x) diff(t.test(x[1:22],x[23:61])$estimate))
  #ciL.mestdagh.est.max[ind.DEt,s]<-apply(normData.max,2,function(x) t.test(x[23:61],x[1:22])$conf.int[1])
  #ciU.mestdagh.est.max[ind.DEt,s]<-apply(normData.max,2,function(x) t.test(x[23:61],x[1:22])$conf.int[2])
  p.mestdagh.est.knn[ind.DEt,s]<-apply(normData.knn,2,function(x) wilcox.test(x[1:22],x[23:61])$p.value)
  
  de.mestdagh.mod.glob[ind.DEt,s]<-apply(normData.mod.glob,2,function(x) diff(t.test(x[1:22],x[23:61])$estimate))
  p.mestdagh.mod.glob[ind.DEt,s]<-apply(normData.mod.glob,2,function(x) wilcox.test(x[1:22],x[23:61])$p.value)
  
  de.mestdagh.mod.glob.ct[ind.DEt,s]<-apply(normData.mod.glob.ct,2,function(x) diff(t.test(x[1:22],x[23:61])$estimate))
  p.mestdagh.mod.glob.ct[ind.DEt,s]<-apply(normData.mod.glob.ct,2,function(x) wilcox.test(x[1:22],x[23:61])$p.value)

  print("Na ttest\n")
  
  cq.during.sim.varL<-replace(cq.during.sim.varL,is.na(cq.during.sim.varL),vLOD)
  # Weighted or not??
  tobit.fit<-tobit(cq.during.sim.varL~mm.unified-1,left=-Inf,right=vLOD)
  tobit.fit.res<-lm(residuals(tobit.fit)^2~miRNA)
  var.tobit.fitted<-(fitted(tobit.fit.res)^{-1})
  tobit.subj.unified.s.weighted<-summary(tobit(cq.during.sim.varL~mm.unified-1,left=-Inf,right=vLOD,weights=var.tobit.fitted))
  tobit.subj.unified.s.varL<-summary(tobit(cq.during.sim.varL~mm.unified-1,left=-Inf,right=vLOD))
  de.tobit.est.varL[ind.DE,s]<-coef(tobit.subj.unified.s.varL)[(nMi+nS):(2*nMi+nS-2),1]
  de.tobit.est.weighted[ind.DE,s]<-coef(tobit.subj.unified.s.weighted)[(nMi+nS):(2*nMi+nS-2),1]
  #std.tobit.est.varL[ind.DE,s]<-coef(tobit.subj.unified.s.varL)[(nMi+nS):(2*nMi+nS-2),2]
  
  p.tobit.est.varL[ind.DE,s]<-coef(tobit.subj.unified.s.varL)[(nMi+nS):(2*nMi+nS-2),4]
  p.tobit.est.weighted[ind.DE,s]<-coef(tobit.subj.unified.s.weighted)[(nMi+nS):(2*nMi+nS-2),4]
  mse.de.tobit.varL[ind.DE,s]<-((de.tobit.est.varL[ind.DE,s]-de.tobit.orig[ind.DE])) ### update 21/02/2012 left squared out
  mse.de.mestdagh.varL[ind.DE,s]<-((de.mestdagh.est.varL[ind.DE,s]-de.mestdagh.orig[ind.DE]))
  mse.de.mestdagh.max[ind.DE,s]<-((de.mestdagh.est.max[ind.DE,s]-de.mestdagh.orig[ind.DE]))
  mse.de.mestdagh.max1[ind.DE,s]<-((de.mestdagh.est.max1[ind.DE,s]-de.mestdagh.orig[ind.DE]))
  mse.de.tobit.weighted[ind.DE,s]<-((de.tobit.est.weighted[ind.DE,s]-de.tobit.orig[ind.DE]))
  mse.de.mestdagh.knn[ind.DE,s]<-((de.mestdagh.est.knn[ind.DE,s]-de.mestdagh.orig[ind.DE]))
  mse.de.mestdagh.mod.glob[ind.DE,s]<-((de.mestdagh.mod.glob[ind.DE,s]-de.mestdagh.orig[ind.DE]))
  mse.de.mestdagh.mod.glob.ct[ind.DE,s]<-((de.mestdagh.mod.glob.ct[ind.DE,s]-de.mestdagh.orig[ind.DE]))
  
  conc.tobit[ind.DE,s]<-ifelse(p.tobit.est.varL[ind.DE,s]<alpha,1,0)
  conc.tobit.weighted[ind.DE,s]<-ifelse(p.tobit.est.weighted[ind.DE,s]<alpha,1,0)
  conc.mestdagh.varL[ind.DE,s]<-ifelse(p.mestdagh.est.varL[ind.DE,s]<alpha,1,0)
  conc.mestdagh.max[ind.DE,s]<-ifelse(p.mestdagh.est.max[ind.DE,s]<alpha,1,0)
  conc.mestdagh.max1[ind.DE,s]<-ifelse(p.mestdagh.est.max1[ind.DE,s]<alpha,1,0)
  conc.mestdagh.knn[ind.DE,s]<-ifelse(p.mestdagh.est.knn[ind.DE,s]<alpha,1,0)
  conc.mestdagh.mod.glob[ind.DE,s]<-ifelse(p.mestdagh.mod.glob[ind.DE,s]<alpha,1,0)
  conc.mestdagh.mod.glob.ct[ind.DE,s]<-ifelse(p.mestdagh.mod.glob.ct[ind.DE,s]<alpha,1,0)
  
  for (i in 2:nS) {
    tobit.beta.est.varL[i-1,s]<-coef(tobit.subj.unified.s.varL)[1,1]+coef(tobit.subj.unified.s.varL)[i+nMi-1,1]
  }
  if (max(tobit.beta.est)>100 | max(tobit.beta.est.varL)>100) {
    print(s)
  }
}
sink()
save.image(file="Wilcox50.RData")