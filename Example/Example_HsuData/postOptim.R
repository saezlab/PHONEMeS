library(igraph)
library(BioNet)
library(PHONEMeS)
setwd("~/Desktop/CamMS_draft/Hsu2011/Cluster/rapa")
load("data4cluster_9.RData")
load("p9_imported.RData")
#This is the folder where all the results will be written
resFolder<-"~/Desktop/CamMS_draft/Hsu2011/Cluster/rapa/"
#These are used to produce the node attribute files and to produce the networks
nodesOnOff<-nodesData(data.On=data.P, dataGMM=dataGMM, pknList=pknList)
#####################
#Write the node attributes
#The node attributes and _NA and data table _DA files are identical for each
#independent round of optimisation so this can be run once
#if multiple independent optimisations are performed 
#(not the case for the edge attributes _EA and network files)
drugsD<-c("rapa - ins")
nodeAttr(targets.On=targets.P, nodesOF=nodesOnOff, drugsD=drugsD, dataGMM=dataGMM, optParam=optParam, resFolder=resFolder)
#####################    
##This produces the edge attributes file
EA<-data.frame(SID=interactions(pknList)$SID, ntag=interactions(pknList)$ntag, f50=opres$FE[match(interactions(pknList)$SID, rownames(opres$FE)),dim(opres$FE)[2]])
write.table(EA, file=paste(resFolder,"EA_p",resN(optParam),".txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
##This writes a maximum input network based on the bin numbers for each edge
iN<-opres$intgNone[, dim(opres$intgNone)[2]]
names(iN)<-rownames(opres$intgNone)

complete.I.mod<-interactions(pknList)
complete.I.mod$ntag<-opres$G1.freq[,dim(opres$G1.freq)[2]]

mI.ntag<-mInw(nwTable=complete.I.mod, targets.On=targets.P, intgNone=iN, nodesOF=nodesOnOff, tol=0)
write.table(mI.ntag, file=paste(resFolder,"ntagMaxIn_PR_p",resN(optParam),".txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE)
##This writes a maximum input network based on the frequency of each edge in the optimised population
complete.I.mod<-interactions(pknList)
complete.I.mod$ntag<-opres$FE[match(interactions(pknList)$SID, rownames(opres$FE)),dim(opres$FE)[2]]
mI.freq<-mInw(nwTable=complete.I.mod, intgNone=iN, nodesOF=nodesOnOff, tol=0, targets.On=targets.P)
write.table(mI.freq, file=paste(resFolder,"maxIn_freq_PR_p",resN(optParam),".txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE)
####################
save.image(paste(resFolder,"p9_processed.RData", sep=""))
###################
#this write in a data file the objects that are needed when multiple independent
#optimisations are performed - ignore if only one optimisation is done
FE.9<-opres$FE
Msize.9<-opres$Msize
BMGsize.9<-opres$BMGsize
BM.s.9<-opres$BM.s
sM.9<-opres$sM
sM.avg.9<-opres$sM.avg
sAll.9<-opres$sAll
G1.freq.9<-opres$G1.freq
G1.flipP.9<-opres[[12]]
G1.AndBinF.9<-opres[[13]]
G1.NoneBinF.9<-opres[[14]]
intgAnd.9<-opres$intgAnd
intgNone.9<-opres$intgNone
save(file=paste(resFolder,"objects_p9.RData", sep=""),
     list=c("intgAnd.9", "intgNone.9","G1.AndBinF.9","G1.flipP.9","FE.9","Msize.9","BMGsize.9","BM.s.9","sM.9","sM.avg.9","sAll.9","G1.freq.9", "G1.NoneBinF.9"))
######################
rm(list=ls())
load("p10_imported.RData")
#This is the folder where all the results will be written
resFolder<-"~/Desktop/CamMS_draft/Hsu2011/Cluster/rapa/"
#this write in a data file the objects that are needed when multiple independent
#optimisations are performed - ignore if only one optimisation is done
FE.10<-opres$FE
Msize.10<-opres$Msize
BMGsize.10<-opres$BMGsize
BM.s.10<-opres$BM.s
sM.10<-opres$sM
sM.avg.10<-opres$sM.avg
sAll.10<-opres$sAll
G1.freq.10<-opres$G1.freq
G1.flipP.10<-opres[[12]]
G1.AndBinF.10<-opres[[13]]
G1.NoneBinF.10<-opres[[14]]
intgAnd.10<-opres$intgAnd
intgNone.10<-opres$intgNone
save(file=paste(resFolder,"objects_p10.RData", sep=""),
     list=c("intgAnd.10", "intgNone.10","G1.AndBinF.10","G1.flipP.10","FE.10","Msize.10","BMGsize.10","BM.s.10","sM.10","sM.avg.10","sAll.10","G1.freq.10", "G1.NoneBinF.10"))
######################
rm(list=ls())
load("p11_imported.RData")
#This is the folder where all the results will be written
resFolder<-"~/Desktop/CamMS_draft/Hsu2011/Cluster/rapa/"
#this write in a data file the objects that are needed when multiple independent
#optimisations are performed - ignore if only one optimisation is done
FE.11<-opres$FE
Msize.11<-opres$Msize
BMGsize.11<-opres$BMGsize
BM.s.11<-opres$BM.s
sM.11<-opres$sM
sM.avg.11<-opres$sM.avg
sAll.11<-opres$sAll
G1.freq.11<-opres$G1.freq
G1.flipP.11<-opres[[12]]
G1.AndBinF.11<-opres[[13]]
G1.NoneBinF.11<-opres[[14]]
intgAnd.11<-opres$intgAnd
intgNone.11<-opres$intgNone
save(file=paste(resFolder,"objects_p11.RData", sep=""),
     list=c("intgAnd.11", "intgNone.11","G1.AndBinF.11","G1.flipP.11","FE.11","Msize.11","BMGsize.11","BM.s.11","sM.11","sM.avg.11","sAll.11","G1.freq.11", "G1.NoneBinF.11"))
######################
rm(list=ls())
load("p12_imported.RData")
#This is the folder where all the results will be written
resFolder<-"~/Desktop/CamMS_draft/Hsu2011/Cluster/rapa/"
#this write in a data file the objects that are needed when multiple independent
#optimisations are performed - ignore if only one optimisation is done
FE.12<-opres$FE
Msize.12<-opres$Msize
BMGsize.12<-opres$BMGsize
BM.s.12<-opres$BM.s
sM.12<-opres$sM
sM.avg.12<-opres$sM.avg
sAll.12<-opres$sAll
G1.freq.12<-opres$G1.freq
G1.flipP.12<-opres[[12]]
G1.AndBinF.12<-opres[[13]]
G1.NoneBinF.12<-opres[[14]]
intgAnd.12<-opres$intgAnd
intgNone.12<-opres$intgNone
save(file=paste(resFolder,"objects_p12.RData", sep=""),
     list=c("intgAnd.12", "intgNone.12","G1.AndBinF.12","G1.flipP.12","FE.12","Msize.12","BMGsize.12","BM.s.12","sM.12","sM.avg.12","sAll.12","G1.freq.12", "G1.NoneBinF.12"))
######################
rm(list=ls())
load("p13_imported.RData")
#This is the folder where all the results will be written
resFolder<-"~/Desktop/CamMS_draft/Hsu2011/Cluster/rapa/"
#this write in a data file the objects that are needed when multiple independent
#optimisations are performed - ignore if only one optimisation is done
FE.13<-opres$FE
Msize.13<-opres$Msize
BMGsize.13<-opres$BMGsize
BM.s.13<-opres$BM.s
sM.13<-opres$sM
sM.avg.13<-opres$sM.avg
sAll.13<-opres$sAll
G1.freq.13<-opres$G1.freq
G1.flipP.13<-opres[[12]]
G1.AndBinF.13<-opres[[13]]
G1.NoneBinF.13<-opres[[14]]
intgAnd.13<-opres$intgAnd
intgNone.13<-opres$intgNone
save(file=paste(resFolder,"objects_p13.RData", sep=""),
     list=c("intgAnd.13", "intgNone.13","G1.AndBinF.13","G1.flipP.13","FE.13","Msize.13","BMGsize.13","BM.s.13","sM.13","sM.avg.13","sAll.13","G1.freq.13", "G1.NoneBinF.13"))
######################
rm(list=ls())
load("p14_imported.RData")
#This is the folder where all the results will be written
resFolder<-"~/Desktop/CamMS_draft/Hsu2011/Cluster/rapa/"
#this write in a data file the objects that are needed when multiple independent
#optimisations are performed - ignore if only one optimisation is done
FE.14<-opres$FE
Msize.14<-opres$Msize
BMGsize.14<-opres$BMGsize
BM.s.14<-opres$BM.s
sM.14<-opres$sM
sM.avg.14<-opres$sM.avg
sAll.14<-opres$sAll
G1.freq.14<-opres$G1.freq
G1.flipP.14<-opres[[12]]
G1.AndBinF.14<-opres[[13]]
G1.NoneBinF.14<-opres[[14]]
intgAnd.14<-opres$intgAnd
intgNone.14<-opres$intgNone
save(file=paste(resFolder,"objects_p14.RData", sep=""),
     list=c("intgAnd.14", "intgNone.14","G1.AndBinF.14","G1.flipP.14","FE.14","Msize.14","BMGsize.14","BM.s.14","sM.14","sM.avg.14","sAll.14","G1.freq.14", "G1.NoneBinF.14"))
