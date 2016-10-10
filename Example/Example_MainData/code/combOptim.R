resFolder<-"~/Desktop/Res"
load(paste(resFolder,"objects_p3.RData", sep=""))
load(paste(resFolder,"objects_p4.RData", sep=""))
load(paste(resFolder,"objects_p5.RData", sep=""))
load(paste(resFolder,"objects_p6.RData", sep=""))
load(paste(resFolder,"objects_p7.RData", sep=""))
load(paste(resFolder,"objects_p8.RData", sep=""))
library(PHONEMeS)
library(BioNet)
library(igraph)
#
#change the number here depending on the number of otimisations that you are combining
colVec<-rainbow(6)
#this number is the final number of generations
fGen<-50
#this is the number of models in each generation
nM<-5000
#matrices of population average model size and best model (of the generation) size
avgMsize<-cbind(rowMeans(Msize.3), 
	rowMeans(Msize.4, na.rm=TRUE), 
	rowMeans(Msize.5, na.rm=TRUE), 
	rowMeans(Msize.6, na.rm=TRUE), 
	rowMeans(Msize.7, na.rm=TRUE), 
	rowMeans(Msize.8, na.rm=TRUE))
BMGsize<-cbind(BMGsize.3, BMGsize.4, BMGsize.5, BMGsize.6, BMGsize.7, BMGsize.8)
#matrices of total score average and best model total score
avgsM_t<-cbind(rowSums(sM.avg.3),
	rowSums(sM.avg.4),
	rowSums(sM.avg.5),
	rowSums(sM.avg.6),
	rowSums(sM.avg.7),
	rowSums(sM.avg.8))
BMGsM_t<-cbind(rowSums(sM.3),
	rowSums(sM.4),
	rowSums(sM.5),
	rowSums(sM.6),
	rowSums(sM.7),
	rowSums(sM.8))
#list of matrices of score average and best model, by drug (1 drug=1 element of the list)
avgsM<-list(drug1=cbind(sM.avg.3[,1],
	sM.avg.4[,1],
	sM.avg.5[,1],
	sM.avg.6[,1],
	sM.avg.7[,1],
	sM.avg.8[,1]))
BMGsM<-list(drug1=cbind(sM.3[,1],
	sM.4[,1],
	sM.5[,1],
	sM.6[,1],
	sM.7[,1],
	sM.8[,1]))
if(!is.null(dim(sM.3)) && dim(sM.3)[2] > 1){	
	for(i in 2:dim(sM.3)[2]){
		avgsM[[i]]<-cbind(sM.avg.3[,i],
			sM.avg.4[,i],
			sM.avg.5[,i],
			sM.avg.6[,i],
			sM.avg.7[,i],
			sM.avg.8[,i])
		BMGsM[[i]]<-cbind(sM.3[,i],
			sM.4[,i],
			sM.5[,i],
			sM.6[,i],
			sM.7[,i],
			sM.8[,i])
	}
	names(avgsM)<-paste(rep("drug", dim(sM.3)[2]), 1:dim(sM.3)[2], sep="")
	names(BMGsM)<-paste(rep("drug", dim(sM.3)[2]), 1:dim(sM.3)[2], sep="")
}
#######Make the plots
pdf(paste(resFolder, "combOptim_plots.pdf", sep=""))
#size plot
plot(avgMsize[,1], type="l", col=colVec[1], main="Size", ylim=c(min(BMGsize),max(BMGsize)))
lines(BMGsize[,1], lty="dashed", col=colVec[1])
for(i in 2:dim(BMGsize)[2]){
	lines(avgMsize[,i], col=colVec[i])
	lines(BMGsize[,i], lty="dashed", col=colVec[i])
}
#total scores plot
plot(BMGsM_t[,1], type="l", ylim=c(min(BMGsM_t, na.rm=TRUE), max(avgsM_t, na.rm=TRUE)), 
	col=colVec[1], xlim=c(0,fGen), 
     xlab="Generation", ylab="Score", lty="dashed", main="Total score")
lines(avgsM_t[,1], col=colVec[1])
for(i in 2:dim(avgsM_t)[2]){
	lines(BMGsM_t[,i], col=colVec[i], lty="dashed")
	lines(avgsM_t[,i], col=colVec[i])
}
#individual drug scores plots
for(l in 1:length(avgsM)){
	plot(BMGsM[[l]][,1], type="l", col=colVec[1], xlim=c(0,fGen),
		ylim=c(min(BMGsM[[l]], na.rm=TRUE), max(avgsM[[l]], na.rm=TRUE)),  
		xlab="Generation", ylab="Score", lty="dashed",
		main=paste(names(avgsM)[l], "scores", sep=" "))
	lines(avgsM[[l]][,1], col=colVec[1])
	for(i in 2:dim(BMGsM[[l]])[2]){
		lines(BMGsM[[l]][,i], col=colVec[i], lty="dashed")
		lines(avgsM[[l]][,i], col=colVec[i])
	}
}
dev.off()
########
####produce the combined networks
#these produce, for each edge, the averaged (across independent optimisations) final 
#bin numbers and frequency in the population
load(paste(resFolder,"data4cluster_3.RData", sep=""))
G1.avg<-rowMeans(cbind(G1.freq.3[,fGen], 
	G1.freq.4[,fGen],
	G1.freq.5[,fGen], 
	G1.freq.6[,fGen], 
	G1.freq.7[,fGen], 
	G1.freq.8[,fGen]))
FE.avg<-rowMeans(cbind(
	FE.3[,fGen], 
	FE.4[,fGen],
	FE.5[,fGen], 
	FE.6[,fGen], 
	FE.7[,fGen], 
	FE.8[,fGen]))
names(G1.avg)<-interactions(pknList)$SID
names(FE.avg)<-rownames(FE.3)
#these are the And and none bin sizes for integrators, and the And flip probabilities
#for intermediates	
intgAnd.avg<-rowMeans(cbind(intgAnd.3[,fGen], 
	intgAnd.4[,fGen],
	intgAnd.5[,fGen], 
	intgAnd.7[,fGen], 
	intgAnd.6[,fGen], 
	intgAnd.8[,fGen]))
intgNone.avg<-rowMeans(cbind(intgNone.3[,fGen], 
	intgNone.4[,fGen],
	intgNone.5[,fGen], 
	intgNone.7[,fGen], 
	intgNone.6[,fGen], 
	intgNone.8[,fGen]))
G1.flipP.avg<-rowMeans(cbind(G1.flipP.3[,fGen], 
	G1.flipP.4[,fGen],
	G1.flipP.5[,fGen], 
	G1.flipP.7[,fGen], 
	G1.flipP.6[,fGen], 
	G1.flipP.8[,fGen]))
#make the elements that tell which nodes to connect in solution networks
nodesOnOff<-nodesData(data.On=data.P, dataGMM=dataGMM, pknList=pknList)
#make the edge attribute file
EA<-data.frame(SID=interactions(pknList)$SID, ntag=G1.avg, f50=FE.avg[interactions(pknList)$SID])
EA$ntag.3<-G1.freq.3[,fGen]
EA$ntag.4<-G1.freq.4[,fGen]
EA$ntag.5<-G1.freq.5[,fGen]
EA$ntag.7<-G1.freq.7[,fGen]
EA$ntag.6<-G1.freq.6[,fGen]
EA$ntag.8<-G1.freq.8[,fGen]
EA$Ffinal.3<-FE.3[match(interactions(pknList)$SID, rownames(FE.3)),fGen]
EA$Ffinal.4<-FE.4[match(interactions(pknList)$SID, rownames(FE.4)),fGen]
EA$Ffinal.5<-FE.5[match(interactions(pknList)$SID, rownames(FE.5)),fGen]
EA$Ffinal.6<-FE.6[match(interactions(pknList)$SID, rownames(FE.6)),fGen]
EA$Ffinal.7<-FE.7[match(interactions(pknList)$SID, rownames(FE.7)),fGen]
EA$Ffinal.8<-FE.8[match(interactions(pknList)$SID, rownames(FE.8)),fGen]
write.table(EA, file=paste(resFolder,"combOptim_EA.txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
###make the max input network by bin numbers
interactions(pknList)$ntag<-G1.avg
iN<-rowMeans(cbind(G1.NoneBinF.3[,fGen], 
	G1.NoneBinF.4[,fGen], 
	G1.NoneBinF.5[,fGen], 
	G1.NoneBinF.6[,fGen], 
	G1.NoneBinF.7[,fGen], 
	G1.NoneBinF.8[,fGen]))
names(iN)<-rownames(G1.NoneBinF.3)
mI.ntag<-mInw(nwTable=interactions(pknList), intgNone=iN, nodesOF=nodesOnOff, tol=0, targets.On=targets.P)
write.table(mI.ntag, file=paste(resFolder,"ntagMaxIn_comb.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
###make the max input networks by frequencies in optimised population - with various tolerances
complete.I.mod<-interactions(pknList)
complete.I.mod$ntag<-FE.avg[interactions(pknList)$SID]
##maxin freq - 0%
mIf.f0<-mInw(nwTable=complete.I.mod, intgNone=intgNone.avg, nodesOF=nodesOnOff, tol=0, targets.On=targets.P)
dim(mIf.f0)#39
write.table(mIf.f0, file=paste(resFolder,"maxInFreq_tol0pc_comb.txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE)
##maxin freq - 1%
mIf.f1<-mInw(nwTable=complete.I.mod, intgNone=intgNone.avg, nodesOF=nodesOnOff, tol=0.01, targets.On=targets.P)
dim(mIf.f1)#41
write.table(mIf.f1, file=paste(resFolder,"maxInFreq_tol1pc_comb.txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE)
##maxin freq - 5%
mIf.f5<-mInw(nwTable=complete.I.mod, intgNone=intgNone.avg, nodesOF=nodesOnOff, tol=0.05, targets.On=targets.P)
dim(mIf.f5)#45
write.table(mIf.f5, file=paste(resFolder,"maxInFreq_tol5pc_comb.txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE)
##maxin freq - 10%
mIf.f10<-mInw(nwTable=complete.I.mod, intgNone=intgNone.avg, nodesOF=nodesOnOff, tol=0.1, targets.On=targets.P)
dim(mIf.f10)#47
write.table(mIf.f10, file=paste(resFolder,"maxInFreq_tol10pc_comb.txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE)
##maxin freq - 20%
mIf.f20<-mInw(nwTable=complete.I.mod, intgNone=intgNone.avg, nodesOF=nodesOnOff, tol=0.2, targets.On=targets.P)
dim(mIf.f20)#58
write.table(mIf.f20, file=paste(resFolder,"maxInFreq_tol20pc_comb.txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE)
####make the max scoring path to single drug targets (based on frequencies)
#this is only an example on one drug and one drug target, it can be repeated for each 
#drug and targets (it finds the maximum frequency path between each node in a set of
#perturbed nodes and a single kinase of your choice
MTORpaths<-pathsD(nwTable=complete.I.mod, nM=nG1(optParam)*nScripts(optParam),nodes2link=unlist(nodesOnOff$On),
drug="MTOR_HUMAN")
write.table(MTORpaths, file=paste(resFolder,"pathsMTOR_freqIcomb.txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE)
#########################
save.image(paste(resFolder,"Res_combined.RData", sep=""))
#########################
nodesP<-rep(NA, length(unique(unlist(nodesOnOff))))
names(nodesP)<-unique(unlist(nodesOnOff))
for(i in 1:length(nodesP)){
  if(names(nodesP)[i] %in% unlist(nodesOnOff$On) && 
       names(nodesP)[i] %in% unlist(nodesOnOff$Off)){
        nodesP[i]<-"B"
  }else{
    if(names(nodesP)[i] %in% unlist(nodesOnOff$On)) nodesP[i]<-"P"
    if(names(nodesP)[i] %in% unlist(nodesOnOff$Off)) nodesP[i]<-"C"
  }
}
write.table(nodesP, file=paste(resFolder,"nodesP_manual.txt", sep=""),
            sep="\t", quote=FALSE, row.names=TRUE)
