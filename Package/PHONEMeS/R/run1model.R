run1model<-function(pknList=pknList, 
           targetsOn=targets.P,
           dataOn=data.P,
           FlipP=rep(0.5, length(intermediates(pknList))),
           AndBinF=rep(0, length(integrators(pknList))),
           NoneBinF=rep(0, length(integrators(pknList))),
           dataGMM=dataGMM,
           optParam){  
    dataOn<-dataBC(dataOn)         
    complete.I<-interactions(pknList)
    sinks<-sinks(pknList)
    nodes<-species(pknList)
    complete.I.start<-interactionsD(pknList)
    intermediates<-intermediates(pknList)
    integrators<-integrators(pknList)
    GMM.res.ID<-IDmap(dataGMM)
    #run the intermediates AND sample
    intermediates.AND<-(runif(length(intermediates), 0, 1) < FlipP)
    if(intgAsintm(optParam) == TRUE){
      integrators.AND<-(runif(length(integrators), 0, 1) < AndBinF)
    }
    #create the model
    model<-createModel(pknList,
                       AndBinF=AndBinF,
                       NoneBinF=NoneBinF,
                       optParam=optParam)                  
    
    scoresT<-0
    scores.dataOn<-dataOn
    intm.inN<-vector("list", length(intermediates))
    for(i in 1:length(intermediates)){
      intm.inN[[i]]<-model[which(model[,"S.cc"] == intermediates[i]),"K.ID"]
    }    
    intg.inN<-vector("list", length(integrators))
    for(i in 1:length(integrators)){
      intg.inN[[i]]<-model[which(model[,"S.cc"] == integrators[i]),"K.ID"]
    }	
    sink.inN<-vector("list", length(sinks))
    for(i in 1:length(sinks)){
      sink.inN[[i]]<-model[which(model[,"S.cc"] == sinks[i]),"K.ID"]
    }
    names(intm.inN)<-intermediates
    names(intg.inN)<-integrators
    names(sink.inN)<-sinks
    intm.inN.l<-unlist(lapply(intm.inN, length))
    sink.inN.l<-unlist(lapply(sink.inN, length))
    intg.inN.l<-unlist(lapply(intg.inN, length))
    intm.inN<-intm.inN[which(intm.inN.l > 0)]
    sink.inN<-sink.inN[which(sink.inN.l > 0)]
    intg.inN<-intg.inN[which(intg.inN.l > 0)]
    intermediates.AND.or<-intermediates.AND
    intermediates.AND<-intermediates.AND[which(intm.inN.l > 0)]
    if(intgAsintm(optParam) == TRUE){
      integrators.AND.or<-integrators.AND
      integrators.AND<-integrators.AND[which(intg.inN.l > 0)]
    }else{
      integrators.AND<-NA
      integrators.AND.or<-NA
    }
    for(tOn in 1:length(targetsOn)){
      nodes.t<-rep(0, length(nodes))
      nodes.t<-data.frame(t0=nodes.t, t1=nodes.t, stringsAsFactors=FALSE)
      rownames(nodes.t)<-nodes
      #set the drug node to 1
      nodes.t[targetsOn[[tOn]],1]<-1
      nodes.t<-simulateT1(nodes.t=nodes.t,
                          intermediates.AND=intermediates.AND,
                          integrators.AND=integrators.AND,
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
                            intermediates.AND=intermediates.AND,
                            integrators.AND=integrators.AND,
                            intm.inN=intm.inN,
                            intg.inN=intg.inN,
                            sink.inN=sink.inN)
        nodes.t[targetsOn[[tOn]],2]<-1
        stopcond<-all(nodes.t[,1] == nodes.t[,2])
        if(count > 30) break
        nodes.t[,1]<-nodes.t[,2]
        count<-count+1
      }
      #compute the score for the target  targetsOn[tOn] 
      #based on the data from dataOn[tOn]
      score<-0  
      for(dt in 1:length(dataOn[[tOn]])){
        score.dt <- sum(as.numeric(dataOn[[tOn]][[dt]][intersect(row.names(dataOn[[tOn]][[dt]]), as.character(IDmap(dataGMM)[IDmap(dataGMM)[, "S.cc"] %in% rownames(nodes.t[which(nodes.t[, 2] == 1), ]), "dataID"])), 1]), na.rm = TRUE)
        hits <- rownames(dataOn[[tOn]][[dt]])[intersect(which(dataOn[[tOn]][[dt]][, 2] == "P"), which(dataOn[[tOn]][[dt]][, 4] == "OK"))]
        dataMhits<-rownames(dataOn[[tOn]][[dt]])[which(rownames(dataOn[[tOn]][[dt]]) %in% as.character(IDmap(dataGMM)[as.character(IDmap(dataGMM)[,"S.cc"]) %in% rownames(nodes.t[which(nodes.t[,2] == 1),]), "dataID"]))]
        hits<-hits[!(hits %in% dataMhits)] 
        score.dt<-score.dt+sum(as.numeric(dataOn[[tOn]][[dt]][hits,1])*-1, na.rm=TRUE)
        scores.dataOn[[tOn]][[dt]]<-score.dt
        score<-score+score.dt
      }
      #now compute the total score
      scoresT<-scoresT+score
    }
    ModelList<-new("ModelList", 
                   model=model, 
                   scores=scoresT, 
                   scoresList=scores.dataOn, 
                   intAnd=intermediates.AND.or, 
                   intgAnd=integrators.AND.or)
    return(ModelList) 
  }
