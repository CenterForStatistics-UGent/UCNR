### First obtain simulation results using the Robustness study scripts
### After these results have been obtained, various figures may be generated by using the code below

source("scripts/auxfunctions.R")

par(mar=c(5.1,4.6,4.1,2.1))
plot.bias.all(idx=1000,KNN=T,MNV1=T,LOD=T,maxY=0.17,-0.1)

par(mar=c(5.1,4.8,4.1,2.1))
plot.var.all(idx=1000,KNN=T,MNV1=T,LOD=T,minY=0,maxY=0.53)

plot.pvalues(miR=7,idx=353,KNN=T,LOD=T,MNV1=T,minY=0,maxP=0.05,maxY=4)

