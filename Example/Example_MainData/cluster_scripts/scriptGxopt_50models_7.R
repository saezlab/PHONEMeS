Mindex<-as.numeric(commandArgs(trailingOnly=TRUE)[2])
Gindex<-as.numeric(commandArgs(trailingOnly=TRUE)[3])
library(igraph)
library(BioNet)
library(PHONEMeS)
#check that the previous generation has been processed
load("data4cluster_7.RData")
if(Gindex != 1){
  genResults<-paste("Results_",resN(optParam),"/G",(Gindex-1),"output.RData", sep="")  
  load(file=genResults)
}else{
  FlipP<-rep(0.5,length(intermediates(pknList)))
  AndBinF<-rep(1, length(integrators(pknList)))
  NoneBinF<-rep(1, length(integrators(pknList)))
}
#if direct interacs from drugs should only be forced in at G1:
#if(Gindex == 1){
#  cstart(OptParam)<-TRUE
#	}else{
#		cstart(OptParam)<-FALSE
#	}
modelsG1<-vector("list", nG1(optParam))
modelsG1.intAnd<-matrix(NA, nrow=nG1(optParam), ncol=length(intermediates(pknList)))
scoresG1<-rep(NA, nG1(optParam))
scoresG1.list<-vector("list", nG1(optParam))
for(g in 1:nG1(optParam)){
  testModel<-run1model(pknList,
                       targetsOn=targets.P,
                       dataOn=data.P,
                       FlipP=FlipP,
                       AndBinF=AndBinF,
                       NoneBinF=NoneBinF,
                       dataGMM=dataGMM,
                       optParam=optParam)
  scoresG1[g]<-scores(testModel)+sizeP(optParam)*(dim(model(testModel))[1]/dim(interactions(pknList))[1])
  modelsG1[[g]]<-model(testModel)
  modelsG1.intAnd[g,]<-intAnd(testModel)
  scoresG1.list[[g]]<-scoresList(testModel)
}
save(file=paste("Results_",resN(optParam),"/modelAndScore_",Mindex,".RData", sep=""), list=c("scoresG1", "modelsG1", "modelsG1.intAnd", "scoresG1.list"))
