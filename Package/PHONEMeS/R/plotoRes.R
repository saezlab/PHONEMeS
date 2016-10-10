#######
plotoRes<-function(opres, data.On=data.P, optParam=optParam){
	pdf(paste("optim_",resN(optParam),".pdf", sep=""))
	plot(rowMeans(opres$Msize, na.rm=TRUE), type="l", main="Average and BM size", 
		xlab="Generations",ylab="Size (# edges)",
		ylim=c(floor(min(opres$BMGsize)),ceiling(max(opres$Msize))))
	lines(opres$BMGsize, lty="dashed", col="grey")
	colVec<-rainbow(dim(opres$sM)[2])
	par(cex=0.7)
	plot(opres$sM[,1], type="l", ylim=c(min(opres$sM, na.rm=TRUE), max(opres$sM, na.rm=TRUE)), 
     	col=colVec[1], xlim=c(0,dim(opres$sM)[1]), xlab="Generation", ylab="Score", lty="dashed", main="BM score")
	if(dim(opres$sM)[2] > 1){
		for(i in 2:dim(opres$sM)[2]){
			lines(opres$sM[,i], col=colVec[i], lty="dashed")
			}
		}	
	legend(legend=unlist(lapply(dataBC(data.On), names)),col=colVec, fill=colVec, x=0, y=max(opres$sM, na.rm=TRUE))
	#
	plot(opres$sM.avg[,1],  type="l", 
		ylim=c(min(opres$sM.avg, na.rm=TRUE), max(opres$sM.avg, na.rm=TRUE)), 
		col=colVec[1], xlim=c(0,dim(opres$sM)[1]), xlab="Generation", ylab="Score", main="Population", lty="solid")
	if(dim(opres$sM.avg)[2] > 1){	
		for(i in 2:dim(opres$sM.avg)[2]){
			lines(opres$sM.avg[,i], col=colVec[i], lty="solid")
			}
		}	
	legend(legend=unlist(lapply(dataBC(data.On), names)),col=colVec, fill=colVec, x=0, y=max(opres$sM.avg, na.rm=TRUE))	
	#
	plot(rowSums(opres$sM.avg, na.rm=TRUE), type="l",xlim=c(0,dim(opres$sM)[1]), 
		ylim=c(min(opres$sM, na.rm=TRUE), max(rowSums(opres$sM.avg), na.rm=TRUE)), 
		xlab="generation", ylab="score", main="Best model and average score, total (no sizeP)")
	lines(rowSums(opres$sM), lty="dashed")
	#
	plot(rowMeans(opres$sAll, na.rm=TRUE), type="l",xlim=c(0,dim(opres$sAll)[1]), 
		ylim=c(min(opres$sAll, na.rm=TRUE), max(opres$sAll, na.rm=TRUE)), 
		xlab="generation", ylab="score", main="Best model and average score, total (with sizeP)")
	lines(apply(opres$sAll, MARGIN=1, FUN=min), lty="dashed")
	nMG<-nScripts(optParam)*nG1(optParam)
	middle<-intersect(which(opres$FE[,dim(opres$FE)[2]] < 0.6*nMG), which(opres$FE[,dim(opres$FE)[2]] > 0.4*nMG))
	top<-which(opres$FE[,dim(opres$FE)[2]] > 0.6*nMG)
	bottom<-which(opres$FE[,dim(opres$FE)[2]] < 0.4*nMG)
	plot(opres$FE[middle[1],]/nMG, ylim=c(0, 1), type="l", 
     	main=paste('Fraction of models where an edge appears (n=',length(middle),')', sep=""), xlab="Generations", ylab="fraction of models")
	for(i in 2:length(middle)){
  		lines(opres$FE[middle[i],]/nMG)
	}
	plot(opres$FE[top[1],]/nMG, ylim=c(0, 1), type="l",
     	main=paste('Fraction of models where an edge appears (n=',length(top),')', sep=""), xlab="Generations", ylab="fraction of models")
	for(i in 2:length(top)){
  		lines(opres$FE[top[i],]/nMG)
	}
	plot(opres$FE[bottom[1],]/nMG, ylim=c(0, 1), type="l", 
     	main=paste('Fraction of models where an edge appears (n=', length(bottom),')', sep=""), xlab="Generations", ylab="fraction of models")
	for(i in 2:length(bottom)){
  		lines(opres$FE[bottom[i],]/nMG)
	}
	dev.off()	
}
