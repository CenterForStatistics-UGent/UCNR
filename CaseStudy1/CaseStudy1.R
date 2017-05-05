## Case study 1
#####################

require(multtest)
require(survival)
require(e1071)
require(AER)
require(maxLik)

### Read data
#############

path<-"";
srcpath<-paste(path,"data/",sep="");
outpath<-paste(path,"results/",sep="");
scriptpath<-paste(path,"scripts/",sep="");

source(paste(scriptpath,"read in neuroblastoma.R",sep=""))
source(paste(scriptpath,"ucnrfunctions.R",sep=""))

names(NB.ok.stack)[1]<-"subject";
ok.stack<-NB.ok.stack

truncated<-35;


### Data preprocessing
######################
### remove miRNAs with less than 20% detects
all.UV.miRNAs<-names(which(table(NB.stack$miRNA,NB.stack$group,NB.stack$determined.ind)[,1,2]>=17))#15
all.UV.miRNAs<-c(all.UV.miRNAs,names(which(table(NB.stack$miRNA,NB.stack$group,NB.stack$determined.ind)[,2,2]>=31)))#33
all.UV.miRNAs<-unique(all.UV.miRNAs)

idx.all.UV<-c()
for (i in 1:length(all.UV.miRNAs)) {
  idx.all.UV<-c(idx.all.UV,grep(all.UV.miRNAs[i],NB.stack$miRNA))
}
NB.stack.filter<-NB.stack[-idx.all.UV,]
NB.data<-NB.stack.filter

### remove determined.ind column, not needed for UCNR model
NB.data <- NB.data[,-5]

### Data formatting
###################

### input for UCNR model needs following columns:
###### sample: sample identifier
###### target: gene identifier
###### Cq: quantification cycle
###### group: two groups to compare

colnames(NB.data)<-c("sample","target","Cq","group")

group<-factor(NB.data$group)
genes<-factor(NB.data$target)
subject<-factor(NB.data$sample)


### fit UCNR
############
fit <- ucnr.fit(data=NB.data, truncated = truncated)

### process fit object to get relevant estimates (e.g. miR-17-92 cluster)
#########################################################################
id.cluster<-c()
id.cluster[1]<-which(rownames(fit)=="mm.unified.extmiRNAhsa-mir-17-3p:groupMYCNSC")
id.cluster[2]<-which(rownames(fit)=="mm.unified.extmiRNAhsa-mir-17-5p:groupMYCNSC")
id.cluster[3]<-which(rownames(fit)=="mm.unified.extmiRNAhsa-mir-18a:groupMYCNSC")
id.cluster[4]<-which(rownames(fit)=="mm.unified.extmiRNAhsa-mir-18asharp:groupMYCNSC")
id.cluster[5]<-which(rownames(fit)=="mm.unified.extmiRNAhsa-mir-19a:groupMYCNSC")
id.cluster[6]<-which(rownames(fit)=="mm.unified.extmiRNAhsa-mir-19b:groupMYCNSC")
id.cluster[7]<-which(rownames(fit)=="mm.unified.extmiRNAhsa-mir-20a:groupMYCNSC")
id.cluster[8]<-which(rownames(fit)=="mm.unified.extmiRNAhsa-mir-92:groupMYCNSC")

p.adjust(fit[,4], method="BH")[id.cluster]