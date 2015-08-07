Gindex<-as.numeric(commandArgs(trailingOnly=TRUE)[2])
library(igraph)
library(BioNet)
library(PHONEMeS)
load("data4cluster_3.RData")
if(Gindex != 1){
  genResults<-paste("Results_",resN(optParam),"/G",(Gindex-1),"output.RData", sep="")
}	
#now I combine all results
Mindex<-1
load(file=paste("Results_",resN(optParam),"/modelAndScore_",Mindex,".RData", sep=""))
c.scoresG1<-scoresG1
c.modelsG1<-modelsG1
c.modelsG1.intAnd<-modelsG1.intAnd
c.scoresG1.list<-scoresG1.list
for(Mindex in 2:nScripts(optParam)){
  load(file=paste("Results_",resN(optParam),"/modelAndScore_",Mindex,".RData", sep=""))
  c.scoresG1<-c(c.scoresG1, scoresG1)
  c.modelsG1<-c(c.modelsG1,modelsG1)
  c.modelsG1.intAnd<-rbind(c.modelsG1.intAnd,modelsG1.intAnd)
  c.scoresG1.list<-c(c.scoresG1.list,scoresG1.list)
}
if(Gindex != 1){
  load(genResults)
  FlipP.p<-FlipP
  AndBinF.p<-AndBinF	
  NoneBinF.p<-NoneBinF
  }else{
  	FlipP.p<-rep(0.5,length(intermediates(pknList)))
  	AndBinF.p<-rep(1, length(integrators(pknList)))
  	NoneBinF.p<-rep(1, length(integrators(pknList)))
	}
###set absTol(optParam)<-FALSE in optParam if you want relative tolerance
if(absTol(optParam) == TRUE){
  models.sel<-c.modelsG1[which(c.scoresG1 < (min(c.scoresG1)+tol(optParam)))]
  models.sel.intAnd<-c.modelsG1.intAnd[which(c.scoresG1 < (min(c.scoresG1)+tol(optParam))),]
	}else{
		models.sel<-c.modelsG1[which(c.scoresG1 < (min(c.scoresG1)+tol(optParam)*abs(min(c.scoresG1))))]
		models.sel.intAnd<-c.modelsG1.intAnd[which(c.scoresG1 < (min(c.scoresG1)+tol(optParam)*abs(min(c.scoresG1)))),]
	}	
if(length(models.sel) == 1){
  models.sel[[2]]<-models.sel[[1]]
  models.sel.intAnd<-rbind(models.sel.intAnd, models.sel.intAnd)
}
models.sel.comb<-models.sel[[1]]
for(i in 2:length(models.sel)){
  models.sel.comb<-rbind(models.sel.comb, models.sel[[i]][!duplicated(models.sel[[i]]),])
}
#####the cap is also a parameter that can be set in optParam
sW<-sinkWeights(pknList=pknList,models.sel.comb=models.sel.comb, models.sel=models.sel, optParam=optParam)
igWeights<-integratorWeights(models.sel=models.sel,
                             pknList=pknList,
                             models.sel.intgAnd=NA,
                             models.sel.comb=models.sel.comb,
                             optParam=optParam)
imWeights<-intermediateWeights(pknList=pknList,
                               models.sel.intAnd=models.sel.intAnd,
                               models.sel=models.sel,
                               models.sel.comb=models.sel.comb,
                             optParam=optParam)        
SID.2paste<-c(sW$SID, igWeights$integrators2paste$SID,  imWeights$interm2paste$SID)     
SID.2paste<-table(SID.2paste)
interactions(pknList)[match(names(SID.2paste), interactions(pknList)[,"SID"]),"ntag"]<-interactions(pknList)[match(names(SID.2paste), interactions(pknList)[,"SID"]),"ntag"]+SID.2paste
FlipP<-(imWeights$AndFlipF+FlipP.p)/2
AndBinF<-igWeights$AndBinF+AndBinF.p
NoneBinF<-igWeights$NoneBinF+NoneBinF.p
#
if(intgAsintm(optParam)){
  AndBinF<-(igWeights$AndBinF+AndBinF.p)/2
}
save(file=paste("Results_",resN(optParam),"/G",Gindex,"output.RData", sep=""), list=c("pknList", "FlipP","AndBinF", "NoneBinF"))
save(file=paste("Results_",resN(optParam),"/G",Gindex,"cResults_",resN(optParam),".RData", sep=""),list=c("c.scoresG1","c.modelsG1","c.modelsG1.intAnd","c.scoresG1.list","pknList"))

