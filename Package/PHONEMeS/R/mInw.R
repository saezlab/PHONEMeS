mInw<-function(nwTable, intgNone, nodesOF, tol, targets.On){
	spIn<-unique(nwTable[,"S.cc"])
	tag<-rep(FALSE, dim(nwTable)[1])
	for(i in 1:length(spIn)){
  		iIn<-which(nwTable[,"S.cc"] == spIn[i])
  		iIns<-nwTable[iIn,"ntag"]
  		if(spIn[i] %in% names(intgNone)){
    		iIn<-iIn[which(iIns >= (max(c(iIns, intgNone[spIn[i]]), na.rm=TRUE)-tol*max(c(iIns, intgNone[spIn[i]]), na.rm=TRUE)))]
  			}else{
    			iIn<-iIn[which(iIns >= (max(iIns, na.rm=TRUE)-tol*max(iIns, na.rm=TRUE)))]   
  			}
  		tag[iIn]<-TRUE
	}
	nwTable<-nwTable[tag,]
	if(any(!is.na(nodesOF$On))){
		nw<-makeNetwork(source=nwTable[,"K.ID"], target=nwTable[,"S.cc"], edgemode="directed")
		nw<-igraph.from.graphNEL(nw)
		targetVn<-match(unlist(nodesOF$On), V(nw)$name)
		if(any(is.na(targetVn))){
			targetVn<-targetVn[!is.na(targetVn)]
		}
		paths.to.data<-data.frame(name=V(nw)$name, lengthres=rep(NA, length(V(nw)$name)))
		for(i in 1:length(V(nw)$name)){
  			Paths<-get.all.shortest.paths(g=nw, from=i, to=targetVn, mode="out", weights=NA)
  			paths.to.data[i,"lengthres"]<-length(Paths$res)
		}
		zeros<-as.character(paths.to.data$name[which(paths.to.data$lengthres == 0)])
		nwTable.nopath<-(nwTable$K.ID %in% zeros)+(nwTable$S.cc %in% zeros)
		nwTable<-nwTable[which(nwTable.nopath == 0),]
	}
	
	#this is used to filter the networks so that the edges in the interactionsD networks get 
	#removed if the are not strictly between perturbed nodes and the drug target under which 
	#they are perturbed
	tAnOn<-targets.On
	for(i in 1:length(targets.On)){
		tAnOn[[i]]<-list(targets=targets.On[[i]], nodesP=unique(unlist(nodesOF$Onlist[[i]])))
	}
	tag<-rep(FALSE, dim(nwTable)[1])
	nM<-max(nwTable$ntag, na.rm=TRUE)
	for(i in 1:dim(nwTable)[1]){
  		if(nwTable$ntag[i] == nM){
    		test<-lapply(tAnOn, function(x){return(nwTable$K.ID[i] %in% x[[1]] && nwTable$S.cc[i] %in% x[[2]])})
    		test<-unlist(test)
    		if(any(test) == FALSE){
      		if(nwTable$K.ID[i] %in% unlist(targets.On)){
      		  tag[i]<-TRUE
      		}
    		}
  		}
	}
	nwTable<-nwTable[!tag,]
	return(nwTable)
}
#####
#if ntag: nwTable=c.I.list$complete.I  (of the last generation, which should be the one in the workspace after the loading function
		  #intgNone=opres$G1.noneBinF[,100]
		  #names(intgNone)=rownames(opres$intgNone)
#if freq: nwTable<-c.I.list$complete.I
		#nwTable$ntag<-opres$FE[match(nwTable$SID, rownames(opres$FE)),100]
		#intgNone=opres$intgNone[,100]
		#names(intgNone)=rownames(opres$intgNone)
#if combined freq:nwTable<-c.I.list$complete.I
				  #nwTable$ntag<-FE.avg[match(nwTable$SID, rownames(FE.19))]
				  #intgNone=intgNone.avg
#if combined ntag:c.I.list$complete.I$ntag<-G1.avg
				 #nwTable<-c.I.list$complete.I
				 #intgNone=G1.noneBinF.avg
#if you don't want to remove the paths not going to data nodes, set nodesOF$On=NA
