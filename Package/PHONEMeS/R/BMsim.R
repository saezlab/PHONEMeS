#####BM - make pR and sim
#this works for any model, not only the best one - in theory I should be able to use this for mI networks as well
BMsim<-function(BM, nodes2link, BM.intAnd, pknList=pknList, targetsOn){
#produce the noPathR version
	BM.nw<-makeNetwork(source=BM[,"K.ID"], target=BM[,"S.cc"], edgemode="directed")
	BM.nw<-igraph.from.graphNEL(BM.nw)
	targetVn<-match(nodes2link, V(BM.nw)$name)
	targetVn<-targetVn[!is.na(targetVn)]
	paths.to.data<-data.frame(name=V(BM.nw)$name, lengthres=rep(NA, length(V(BM.nw)$name)))
	for(i in 1:length(V(BM.nw)$name)){
 		Paths<-get.all.shortest.paths(g=BM.nw, from=i, to=targetVn, mode="out", weights=NA)
  		paths.to.data[i,"lengthres"]<-length(Paths$res)
	}
	zeros<-as.character(paths.to.data$name[which(paths.to.data$lengthres == 0)])
	BM.nopath<-(BM$K.ID %in% zeros)+(BM$S.cc %in% zeros)
	BM.nopathR<-BM[which(BM.nopath == 0),]
#now simulate
	model<-BM
	intermediates.t.AND<-BM.intAnd
	intm.inN<-vector("list", length(intermediates(pknList)))
	for(i in 1:length(c.I.list$intermediates.t)){
  		intm.inN[[i]]<-model[which(model[,"S.cc"] == intermediates(pknList)[i]),"K.ID"]
	}    
	intg.inN<-vector("list", length(c.I.list$integrators))
	for(i in 1:length(c.I.list$integrators)){
  		intg.inN[[i]]<-model[which(model[,"S.cc"] == integrators(pknList)[i]),"K.ID"]
	}  
	sink.inN<-vector("list", length(c.I.list$sinks))
	for(i in 1:length(c.I.list$sinks)){
  		sink.inN[[i]]<-model[which(model[,"S.cc"] == sinks(pknList)[i]),"K.ID"]
	}
	names(intm.inN)<-intermediates(pknList)
	names(intg.inN)<-integrators(pknList)
	names(sink.inN)<-sinks(pknList)
	intm.inN.l<-unlist(lapply(intm.inN, length))
	sink.inN.l<-unlist(lapply(sink.inN, length))
	intg.inN.l<-unlist(lapply(intg.inN, length))
	intm.inN<-intm.inN[which(intm.inN.l > 0)]
	sink.inN<-sink.inN[which(sink.inN.l > 0)]
	intg.inN<-intg.inN[which(intg.inN.l > 0)]
	intermediates.t.AND.or<-intermediates.t.AND
	intermediates.t.AND<-intermediates.t.AND[which(intm.inN.l > 0)]
	tOn.res<-vector("list", length(targetsOn))
	for(tOn in 1:length(targetsOn)){
  		nodes.t<-rep(0, length(species(pknList)))
  		nodes.t<-data.frame(t0=nodes.t, t1=nodes.t, stringsAsFactors=FALSE)
  		rownames(nodes.t)<-c.I.list$nodes
  		#set the drug node to 1
  		nodes.t[targetsOn[[tOn]],1]<-1
  		nodes.t<-simulateT1(nodes.t=nodes.t, 
                      intermediates.t.AND=intermediates.t.AND,
                      integrators.t.AND=NA,
                      intm.inN=intm.inN,
                      intg.inN=intg.inN,
                      sink.inN=sink.inN)
  		stopcond<-all(nodes.t[,1] == nodes.t[,2])
  		nodes.t[,1]<-nodes.t[,2]
  		nodes.t[targetsOn[[tOn]],1]<-1
  		count<-1
  		while(!stopcond){
    		print(count)
    		nodes.t<-simulateT1(nodes.t=nodes.t,
                        intermediates.AND=intermediates.t.AND,
                        integrators.AND=NA,
                        intm.inN=intm.inN,
                        intg.inN=intg.inN,
                        sink.inN=sink.inN)
    		nodes.t[targetsOn[[tOn]],2]<-1
    		stopcond<-all(nodes.t[,1] == nodes.t[,2])
    		if(count > 30) break
    			nodes.t[,1]<-nodes.t[,2]
    			count<-count+1
  			}
  		tOn.res[[tOn]]<-nodes.t
	} 
	
	BM.sim<-data.frame(tOn.res[[1]][,2])
	rownames(BM.sim)<-rownames(tOn.res[[1]])
	for(i in 2:length(tOn.res)){
		BM.sim<-cbind(BM.sim, tOn.res[[i]][,2])
	}
	colnames(BM.sim)<-names(targetsOn)
	return(list(BM.nopathR=BM.nopathR, BM.sim=BM.sim))	
}
