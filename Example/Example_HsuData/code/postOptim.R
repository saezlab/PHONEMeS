library(igraph)
library(BioNet)
library(PHONEMeS)

# Set working directory to directory of this script (in RStudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# If you are using this file as Source, use:
# setwd(getSrcDirectory()[1])

#### RAPA ##################

load("../cluster_scripts/rapa/data4cluster_9.RData")
load("../cluster_results/rapa/p9_imported.RData")

#This is the folder where all the results will be written
resFolder<-"../result_network/rapa/"

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
rm(list=setdiff(ls(), "resFolder"))

load("../cluster_results/rapa/p10_imported.RData")

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
rm(list=setdiff(ls(), "resFolder"))

load("../cluster_results/rapa/p11_imported.RData")
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
rm(list=setdiff(ls(), "resFolder"))

load("../cluster_results/rapa/p12_imported.RData")
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
rm(list=setdiff(ls(), "resFolder"))

load("../cluster_results/rapa/p13_imported.RData")
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
rm(list=setdiff(ls(), "resFolder"))

load("../cluster_results/rapa/p14_imported.RData")
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

rm(list=ls())
### TORIN ##############

load("../cluster_scripts/torin/data4cluster_3.RData")
load("../cluster_results/torin/p3_imported.RData")

#This is the folder where all the results will be written
resFolder<-"../result_network/torin/"

#These are used to produce the node attribute files and to produce the networks
nodesOnOff<-nodesData(data.On=data.P, dataGMM=dataGMM, pknList=pknList)
#####################
#Write the node attributes
#The node attributes and _NA and data table _DA files are identical for each
#independent round of optimisation so this can be run once
#if multiple independent optimisations are performed 
#(not the case for the edge attributes _EA and network files)
drugsD<-c("tor - ins")
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
save.image(paste(resFolder,"p3_processed.RData", sep=""))
###################
#this write in a data file the objects that are needed when multiple independent
#optimisations are performed - ignore if only one optimisation is done
FE.3<-opres$FE
Msize.3<-opres$Msize
BMGsize.3<-opres$BMGsize
BM.s.3<-opres$BM.s
sM.3<-opres$sM
sM.avg.3<-opres$sM.avg
sAll.3<-opres$sAll
G1.freq.3<-opres$G1.freq
G1.flipP.3<-opres[[12]]
G1.AndBinF.3<-opres[[13]]
G1.NoneBinF.3<-opres[[14]]
intgAnd.3<-opres$intgAnd
intgNone.3<-opres$intgNone
save(file=paste(resFolder,"objects_p3.RData", sep=""),
     list=c("intgAnd.3", "intgNone.3","G1.AndBinF.3","G1.flipP.3","FE.3","Msize.3","BMGsize.3","BM.s.3","sM.3","sM.avg.3","sAll.3","G1.freq.3", "G1.NoneBinF.3"))
######################
rm(list=setdiff(ls(), "resFolder"))

load("../cluster_results/torin/p4_imported.RData")

FE.4<-opres$FE
Msize.4<-opres$Msize
BMGsize.4<-opres$BMGsize
BM.s.4<-opres$BM.s
sM.4<-opres$sM
sM.avg.4<-opres$sM.avg
sAll.4<-opres$sAll
G1.freq.4<-opres$G1.freq
G1.flipP.4<-opres[[12]]
G1.AndBinF.4<-opres[[13]]
G1.NoneBinF.4<-opres[[14]]
intgAnd.4<-opres$intgAnd
intgNone.4<-opres$intgNone
save(file=paste(resFolder,"objects_p4.RData", sep=""),
     list=c("intgAnd.4", "intgNone.4","G1.AndBinF.4","G1.flipP.4","FE.4","Msize.4","BMGsize.4","BM.s.4","sM.4","sM.avg.4","sAll.4","G1.freq.4", "G1.NoneBinF.4"))
######################
rm(list=setdiff(ls(), "resFolder"))

load("../cluster_results/torin/p5_imported.RData")

FE.5<-opres$FE
Msize.5<-opres$Msize
BMGsize.5<-opres$BMGsize
BM.s.5<-opres$BM.s
sM.5<-opres$sM
sM.avg.5<-opres$sM.avg
sAll.5<-opres$sAll
G1.freq.5<-opres$G1.freq
G1.flipP.5<-opres[[12]]
G1.AndBinF.5<-opres[[13]]
G1.NoneBinF.5<-opres[[14]]
intgAnd.5<-opres$intgAnd
intgNone.5<-opres$intgNone
save(file=paste(resFolder,"objects_p5.RData", sep=""),
     list=c("intgAnd.5", "intgNone.5","G1.AndBinF.5","G1.flipP.5","FE.5","Msize.5","BMGsize.5","BM.s.5","sM.5","sM.avg.5","sAll.5","G1.freq.5", "G1.NoneBinF.5"))
######################
rm(list=setdiff(ls(), "resFolder"))

load("../cluster_results/torin/p6_imported.RData")
FE.6<-opres$FE
Msize.6<-opres$Msize
BMGsize.6<-opres$BMGsize
BM.s.6<-opres$BM.s
sM.6<-opres$sM
sM.avg.6<-opres$sM.avg
sAll.6<-opres$sAll
G1.freq.6<-opres$G1.freq
G1.flipP.6<-opres[[12]]
G1.AndBinF.6<-opres[[13]]
G1.NoneBinF.6<-opres[[14]]
intgAnd.6<-opres$intgAnd
intgNone.6<-opres$intgNone
save(file=paste(resFolder,"objects_p6.RData", sep=""),
     list=c("intgAnd.6", "intgNone.6","G1.AndBinF.6","G1.flipP.6","FE.6","Msize.6","BMGsize.6","BM.s.6","sM.6","sM.avg.6","sAll.6","G1.freq.6", "G1.NoneBinF.6"))
######################
rm(list=setdiff(ls(), "resFolder"))

load("../cluster_results/torin/p7_imported.RData")
FE.7<-opres$FE
Msize.7<-opres$Msize
BMGsize.7<-opres$BMGsize
BM.s.7<-opres$BM.s
sM.7<-opres$sM
sM.avg.7<-opres$sM.avg
sAll.7<-opres$sAll
G1.freq.7<-opres$G1.freq
G1.flipP.7<-opres[[12]]
G1.AndBinF.7<-opres[[13]]
G1.NoneBinF.7<-opres[[14]]
intgAnd.7<-opres$intgAnd
intgNone.7<-opres$intgNone
save(file=paste(resFolder,"objects_p7.RData", sep=""),
     list=c("intgAnd.7", "intgNone.7","G1.AndBinF.7","G1.flipP.7","FE.7","Msize.7","BMGsize.7","BM.s.7","sM.7","sM.avg.7","sAll.7","G1.freq.7", "G1.NoneBinF.7"))
######################
rm(list=setdiff(ls(), "resFolder"))

load("../cluster_results/torin/p8_imported.RData")
FE.8<-opres$FE
Msize.8<-opres$Msize
BMGsize.8<-opres$BMGsize
BM.s.8<-opres$BM.s
sM.8<-opres$sM
sM.avg.8<-opres$sM.avg
sAll.8<-opres$sAll
G1.freq.8<-opres$G1.freq
G1.flipP.8<-opres[[12]]
G1.AndBinF.8<-opres[[13]]
G1.NoneBinF.8<-opres[[14]]
intgAnd.8<-opres$intgAnd
intgNone.8<-opres$intgNone
save(file=paste(resFolder,"objects_p8.RData", sep=""),
     list=c("intgAnd.8", "intgNone.8","G1.AndBinF.8","G1.flipP.8","FE.8","Msize.8","BMGsize.8","BM.s.8","sM.8","sM.avg.8","sAll.8","G1.freq.8", "G1.NoneBinF.8"))
