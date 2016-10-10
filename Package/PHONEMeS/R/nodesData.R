nodesData<-function(data.On, dataGMM, pknList){
	onNodes<-dataBC(data.On)
	offNodes<-dataBC(data.On)
	for(i in 1:length(dataBC(data.On))){
		for(j in 1:length(dataBC(data.On)[[i]])){
			onNodes[[i]][[j]]<-unique(as.character(IDmap(dataGMM)[match(rownames(dataBC(data.On)[[i]][[j]])[intersect(which(dataBC(data.On)[[i]][[j]][,2] == "P"), which(dataBC(data.On)[[i]][[j]][,4] == "OK"))], as.character(IDmap(dataGMM)[,"dataID"])), "S.cc"]))
			offNodes[[i]][[j]]<-unique(as.character(IDmap(dataGMM)[match(rownames(dataBC(data.On)[[i]][[j]])[intersect(which(dataBC(data.On)[[i]][[j]][,2] == "C"), which(dataBC(data.On)[[i]][[j]][,4] == "OK"))], as.character(IDmap(dataGMM)[,"dataID"])), "S.cc"]))
			onNodes[[i]][[j]]<-onNodes[[i]][[j]][onNodes[[i]][[j]] %in% interactions(pknList)[,"S.cc"]]
			offNodes[[i]][[j]]<-offNodes[[i]][[j]][offNodes[[i]][[j]] %in% interactions(pknList)[,"S.cc"]]
		}
	}
	onNodes.2<-lapply(onNodes, function(x){unique(unlist(x))})
	offNodes.2<-lapply(offNodes, function(x){unique(unlist(x))})
	return(list(On=onNodes.2, Off=offNodes.2, Onlist=onNodes, Offlist=offNodes))
}	