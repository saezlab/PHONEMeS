oRes<-function(data.On, optParam, pknList){
	#nG is the number of generations
	#number of models within tolerance
	nM<-rep(NA, nG(optParam))
	#number of drugs (not conditions! a condition can contain >1 drug)
	nDrugs<-sum(unlist(lapply(dataBC(data.On), length)))
	#score of best model of generation, by drug
	sM<-matrix(NA, ncol=nDrugs, nrow=nG(optParam))
	#average score of generation, by drug
	sM.avg<-matrix(NA, ncol=nDrugs, nrow=nG(optParam))
	#total scores of all models, at each generation
	sAll<-matrix(NA, nrow=nG(optParam), ncol=nScripts(optParam)*nG1(optParam))
	#size of all models at each generation
	Msize<-matrix(NA, nrow=nG(optParam), ncol=nScripts(optParam)*nG1(optParam))
	#size of best model
	BMGsize<-rep(NA, nG(optParam))
	#total score of best model, with size penalty
	BM.s<-rep(NA, nG(optParam))
	#freq of each edge in the generation
	FE<-matrix(NA, ncol=nG(optParam), nrow=dim(interactions(pknList))[1])
	rownames(FE)<-interactions(pknList)[,"SID"]
	#none probability at each generation, for all integrator nodes
	intgNone<-matrix(NA, ncol=nG(optParam), nrow=length(integrators(pknList)))
	rownames(intgNone)<-integrators(pknList)
	#and probability at each generation, for all integrator nodes
	intgAnd<-matrix(NA, ncol=nG(optParam), nrow=length(integrators(pknList)))
	rownames(intgAnd)<-integrators(pknList)
	for(i in 1:nG(optParam)){
  		load(paste("Results_",resN(optParam),"/G",i,"cResults_",resN(optParam),".RData", sep=""))
 		BMGsize[i]<-dim(c.modelsG1[[which(c.scoresG1 == min(c.scoresG1))[1]]])[1]
 		Msize[i,]<-unlist(lapply(c.modelsG1, function(x){return(dim(x)[1])}))	
  		if(tol(optParam) < 1){
  			b<-(min(c.scoresG1, na.rm=TRUE)+tol(optParam)*abs(min(c.scoresG1, na.rm=TRUE)))
  		}else{
  			b<-min(c.scoresG1, na.rm=TRUE)+tol(optParam)
  			}
  		nM[i]<-length(which(c.scoresG1 < b))
  		BM.s[i]<-min(c.scoresG1)
  		sM[i,]<-unlist(c.scoresG1.list[[which(c.scoresG1 == min(c.scoresG1))[1]]])
  		tempList<-lapply(c.scoresG1.list, unlist)
  		for(j in 1:nDrugs){
  				sM.avg[i,j]<-mean(unlist(lapply(tempList, function(x){return(x[j])})), na.rm=TRUE)
  			}
  		sAll[i,]<-c.scoresG1
  		sAll[i,]<-c.scoresG1
  		SIDsG<-lapply(c.modelsG1, function(x){return(unique(x$SID))})
  		SIDsG<-unlist(SIDsG)
  		SIDsG<-table(SIDsG)
  		FE[names(SIDsG),i]<-SIDsG  
  		for(j in 1:length(integrators(pknList))){
    		In<-lapply(c.modelsG1, function(x){integrators(pknList)[j] %in% x$S.cc})
    		And<-lapply(c.modelsG1, function(x){return(length(which(x$S.cc == integrators(pknList)[j])) > 1) })
    		intgNone[j,i]<-sum(!unlist(In))
    		intgAnd[j,i]<-sum(unlist(And))
  		}
	}
	load(paste("Results_",resN(optParam),"/G1output.RData", sep=""))
	G1.freq<-interactions(pknList)[,"ntag"]
	G1.flipP<-FlipP
	G1.AndBinF<-AndBinF
	G1.NoneBinF<-NoneBinF
	for(i in 2:nG(optParam)){
  		load(file=paste("Results_", resN(optParam),"/G",i,"output.RData",sep=""))
  		G1.freq<-cbind(G1.freq, interactions(pknList)[,"ntag"])
  		G1.flipP<-cbind(G1.flipP,FlipP)
  		G1.AndBinF<-cbind(G1.AndBinF,AndBinF)
  		G1.NoneBinF<-cbind(G1.NoneBinF,NoneBinF)
	}
	return(list(nM=nM, sM=sM, sM.avg=sM.avg, sAll=sAll, Msize=Msize, BMGsize=BMGsize, 
		BM.s=BM.s,FE=FE,intgNone=intgNone, intgAnd=intgAnd, 
		G1.freq=G1.freq, G1.flipP=G1.flipP, G1.AndBinF=G1.AndBinF, G1.noneBinF=G1.NoneBinF))
}
