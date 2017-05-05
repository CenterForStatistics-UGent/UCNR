#######################
#Generalized Read Code#
#######################

reorganize<-function(data,nam=NULL,addit=NULL) {
	out<-stack(data);
	if (!is.null(addit)) {
		out<-cbind(out,addit);
	}
	if (!is.null(nam)) {
		names(out)<-nam;
	}

	return(out);
}

search.miRNA<-function(nam, nam.List) {
  idx<-c()
  for (i in 1:length(nam)) {
    if (is.na(grep(nam[i],nam.List))==FALSE) {
      idx[i]<-grep(nam[i],nam.List);
    }
  }
  return (idx);
}

##########################
# Read In Neuroblastoma###
##########################

NB<-read.table(file=paste(srcpath,"validation.txt",sep=""),header=T,sep="\t",na.strings="undetermined")
NB<-as.data.frame(NB)
names(NB)<-c("NB","miRNA","Ct")
miRNA.names<-levels(NB$miRNA) #448
NB.subject<-levels(NB$NB)     #61
# control : as.numeric(NB$Ct)

#metadata
samples<-read.csv(file=paste(srcpath,"samples.csv",sep=""),header=F,sep=",",na.strings="undetermined")
samples<-samples[,-c(2,4:18)]
names(samples)<-c("NB","group")

#merge data with metadata
NB.stack<-merge(NB,samples)
NB.stack<-NB.stack[with(NB.stack, order(miRNA)), ]

# K = 35
threshold<-35;
# expressed values
determined.ind<-as.factor(ifelse(NB.stack$Ct>threshold,1,0))
NB.stack<-cbind(NB.stack,determined.ind)

# Data reorganization for SAM
NB.mat<-reshape(NB,idvar="NB",direction="wide",timevar="miRNA")
rownames(NB.mat)<-NB.mat[,1]
NB.mat<-NB.mat[,-1]
colnames(NB.mat)<-miRNA.names

ordered.samples<-samples[with(samples,order(NB)),]
mna<-grep("MNA",ordered.samples$group)
mycnsc<-grep("MYCNSC",ordered.samples$group)


NB.ok.stack<-NB.stack[NB.stack$determined.ind==0,]
NB.undet.stack<-NB.stack[NB.stack$determined.ind==1,]

NB.ok.stack$miRNA<-factor(NB.ok.stack$miRNA)
NB.undet.stack$miRNA<-factor(NB.undet.stack$miRNA)


# Data frame 'NB' : Original Data frame
# Stack 'NB.stack' : Stacked data
# Stack 'NB.ok.stack' : Stacked data with expressed Cq values
# Stack 'NB.undet.stack' : Stacked data with unexpressed Cq values

