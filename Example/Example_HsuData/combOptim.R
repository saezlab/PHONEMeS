resFolder<-"~/Desktop/CamMS_draft/Hsu2011/Cluster/rapa/"
load(paste(resFolder,"objects_p9.RData", sep=""))
load(paste(resFolder,"objects_p10.RData", sep=""))
load(paste(resFolder,"objects_p11.RData", sep=""))
load(paste(resFolder,"objects_p12.RData", sep=""))
load(paste(resFolder,"objects_p13.RData", sep=""))
load(paste(resFolder,"objects_p14.RData", sep=""))
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
avgMsize<-cbind(rowMeans(Msize.9), 
	rowMeans(Msize.10, na.rm=TRUE), 
	rowMeans(Msize.11, na.rm=TRUE), 
	rowMeans(Msize.12, na.rm=TRUE), 
	rowMeans(Msize.13, na.rm=TRUE), 
	rowMeans(Msize.14, na.rm=TRUE))
BMGsize<-cbind(BMGsize.9, BMGsize.10, BMGsize.11, BMGsize.12, BMGsize.13, BMGsize.14)
#matrices of total score average and best model total score
avgsM_t<-cbind(rowSums(sM.avg.9),
	rowSums(sM.avg.10),
	rowSums(sM.avg.11),
	rowSums(sM.avg.12),
	rowSums(sM.avg.13),
	rowSums(sM.avg.14))
BMGsM_t<-cbind(rowSums(sM.9),
	rowSums(sM.10),
	rowSums(sM.11),
	rowSums(sM.12),
	rowSums(sM.13),
	rowSums(sM.14))
#list of matrices of score average and best model, by drug (1 drug=1 element of the list)
avgsM<-list(drug1=cbind(sM.avg.9[,1],
	sM.avg.10[,1],
	sM.avg.11[,1],
	sM.avg.12[,1],
	sM.avg.13[,1],
	sM.avg.14[,1]))
BMGsM<-list(drug1=cbind(sM.9[,1],
	sM.10[,1],
	sM.11[,1],
	sM.12[,1],
	sM.13[,1],
	sM.14[,1]))
if(!is.null(dim(sM.9)) && dim(sM.9)[2] > 1){	
	for(i in 2:dim(sM.9)[2]){
		avgsM[[i]]<-cbind(sM.avg.9[,i],
			sM.avg.10[,i],
			sM.avg.11[,i],
			sM.avg.12[,i],
			sM.avg.13[,i],
			sM.avg.14[,i])
		BMGsM[[i]]<-cbind(sM.9[,i],
			sM.10[,i],
			sM.11[,i],
			sM.12[,i],
			sM.13[,i],
			sM.14[,i])
	}
	names(avgsM)<-paste(rep("drug", dim(sM.9)[2]), 1:dim(sM.9)[2], sep="")
	names(BMGsM)<-paste(rep("drug", dim(sM.9)[2]), 1:dim(sM.9)[2], sep="")
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
load(paste(resFolder,"data4cluster_9.RData", sep=""))
G1.avg<-rowMeans(cbind(G1.freq.9[,fGen], 
	G1.freq.10[,fGen],
	G1.freq.11[,fGen], 
	G1.freq.12[,fGen], 
	G1.freq.13[,fGen], 
	G1.freq.14[,fGen]))
FE.avg<-rowMeans(cbind(
	FE.9[,fGen], 
	FE.10[,fGen],
	FE.11[,fGen], 
	FE.12[,fGen], 
	FE.13[,fGen], 
	FE.14[,fGen]))
names(G1.avg)<-interactions(pknList)$SID
names(FE.avg)<-rownames(FE.9)
#these are the And and none bin sizes for integrators, and the And flip probabilities
#for intermediates	
intgAnd.avg<-rowMeans(cbind(intgAnd.9[,fGen], 
	intgAnd.10[,fGen],
	intgAnd.11[,fGen], 
	intgAnd.13[,fGen], 
	intgAnd.12[,fGen], 
	intgAnd.14[,fGen]))
intgNone.avg<-rowMeans(cbind(intgNone.9[,fGen], 
	intgNone.10[,fGen],
	intgNone.11[,fGen], 
	intgNone.13[,fGen], 
	intgNone.12[,fGen], 
	intgNone.14[,fGen]))
G1.flipP.avg<-rowMeans(cbind(G1.flipP.9[,fGen], 
	G1.flipP.10[,fGen],
	G1.flipP.11[,fGen], 
	G1.flipP.13[,fGen], 
	G1.flipP.12[,fGen], 
	G1.flipP.14[,fGen]))
#make the elements that tell which nodes to connect in solution networks
nodesOnOff<-nodesData(data.On=data.P, dataGMM=dataGMM, pknList=pknList)
#make the edge attribute file
EA<-data.frame(SID=interactions(pknList)$SID, ntag=G1.avg, f50=FE.avg[interactions(pknList)$SID])
EA$ntag.9<-G1.freq.9[,fGen]
EA$ntag.10<-G1.freq.10[,fGen]
EA$ntag.11<-G1.freq.11[,fGen]
EA$ntag.13<-G1.freq.13[,fGen]
EA$ntag.12<-G1.freq.12[,fGen]
EA$ntag.14<-G1.freq.14[,fGen]
EA$Ffinal.9<-FE.9[match(interactions(pknList)$SID, rownames(FE.9)),fGen]
EA$Ffinal.10<-FE.10[match(interactions(pknList)$SID, rownames(FE.10)),fGen]
EA$Ffinal.11<-FE.11[match(interactions(pknList)$SID, rownames(FE.11)),fGen]
EA$Ffinal.12<-FE.12[match(interactions(pknList)$SID, rownames(FE.12)),fGen]
EA$Ffinal.13<-FE.13[match(interactions(pknList)$SID, rownames(FE.13)),fGen]
EA$Ffinal.14<-FE.14[match(interactions(pknList)$SID, rownames(FE.14)),fGen]
write.table(EA, file=paste(resFolder,"combOptim_EA.txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
###make the max input network by bin numbers
interactions(pknList)$ntag<-G1.avg
iN<-rowMeans(cbind(G1.NoneBinF.9[,fGen], 
	G1.NoneBinF.10[,fGen], 
	G1.NoneBinF.11[,fGen], 
	G1.NoneBinF.12[,fGen], 
	G1.NoneBinF.13[,fGen], 
	G1.NoneBinF.14[,fGen]))
names(iN)<-rownames(G1.NoneBinF.9)
mI.ntag<-mInw(nwTable=interactions(pknList), intgNone=iN, nodesOF=nodesOnOff, tol=0, targets.On=targets.P)
write.table(mI.ntag, file=paste(resFolder,"ntagMaxIn_comb.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
###make the max input networks by frequencies in optimised population - with various tolerances
complete.I.mod<-interactions(pknList)
complete.I.mod$ntag<-FE.avg[interactions(pknList)$SID]
##maxin freq - 0%
mIf.f0<-mInw(nwTable=complete.I.mod, intgNone=intgNone.avg, nodesOF=nodesOnOff, tol=0, targets.On=targets.P)
dim(mIf.f0)#34
write.table(mIf.f0, file=paste(resFolder,"maxInFreq_tol0pc_comb.txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE)
##maxin freq - 1%
mIf.f1<-mInw(nwTable=complete.I.mod, intgNone=intgNone.avg, nodesOF=nodesOnOff, tol=0.01, targets.On=targets.P)
dim(mIf.f1)#34
write.table(mIf.f1, file=paste(resFolder,"maxInFreq_tol1pc_comb.txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE)
##maxin freq - 5%
mIf.f5<-mInw(nwTable=complete.I.mod, intgNone=intgNone.avg, nodesOF=nodesOnOff, tol=0.05, targets.On=targets.P)
dim(mIf.f5)#35
write.table(mIf.f5, file=paste(resFolder,"maxInFreq_tol5pc_comb.txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE)
##maxin freq - 10%
mIf.f10<-mInw(nwTable=complete.I.mod, intgNone=intgNone.avg, nodesOF=nodesOnOff, tol=0.1, targets.On=targets.P)
dim(mIf.f10)#37
write.table(mIf.f10, file=paste(resFolder,"maxInFreq_tol10pc_comb.txt", sep=""),
            sep="\t", quote=FALSE, row.names=FALSE)
##maxin freq - 20%
mIf.f20<-mInw(nwTable=complete.I.mod, intgNone=intgNone.avg, nodesOF=nodesOnOff, tol=0.2, targets.On=targets.P)
dim(mIf.f20)#46
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
