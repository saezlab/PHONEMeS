dataBycond <-function(
  dataGMM=dataGMM, bg=bg, scaled=TRUE,
  rowBycond=list(cond1=c("MTOR1 - Control2","MTOR2 - Control2")))
  {
  #create empty object
  dBc<-vector("list", length(rowBycond))
  names(dBc)<-names(rowBycond)
  #look for peptides that map to sites in the network
  dataReachable<-as.character(IDmap(dataGMM)[as.character(IDmap(dataGMM)[,"S.cc"]) %in% species(bg), "dataID"])
  dataReachable<-dataReachable[!is.na(dataReachable)] 
  dataReachable<-as.character(dataReachable) 
  for(c in 1:length(rowBycond)){
    dBc[[c]]<-vector("list", length(rowBycond[[c]]))
    names(dBc[[c]])<-rowBycond[[c]]
    for(dc in 1:length(rowBycond[[c]])){
      dl<-lapply(res(dataGMM), function(x){return(x[rowBycond[[c]][dc],])})
      DM<-matrix(NA, nrow=length(dl), ncol=4)
      for(i in 1:length(dl)){
        DM[i,]<-dl[[i]]
      }
      rownames(DM)<-names(dl)
      if(any(DM[,1] == Inf, na.rm=T)){
        DM[which(DM[,1] == Inf),1]<-max(as.numeric(DM[-which(DM[,1] == Inf),1]), na.rm=TRUE)
      }
      if(any(DM[,1] == -Inf, na.rm=T)){
        DM[which(DM[,1] == -Inf),1]<-min(as.numeric(DM[-which(DM[,1] == -Inf),1]), na.rm=TRUE)
      }
      DM<-DM[(rownames(DM) %in% dataReachable),]
      dBc[[c]][[dc]]<-DM
      colnames(dBc[[c]][[dc]])<-colnames(res(dataGMM)[[1]])
    }
  }
  if(scaled == TRUE){
    for(i in 1:length(dBc)){
      for(j in 1:length(dBc[[i]])){
        maxP<-max(as.numeric(dBc[[i]][[j]][which(!is.na(dBc[[i]][[j]][,4])),1]), na.rm=TRUE)
        maxN<-abs(min(as.numeric(dBc[[i]][[j]][which(!is.na(dBc[[i]][[j]][,4])),1]), na.rm=TRUE))
        dBc[[i]][[j]][which(as.numeric(dBc[[i]][[j]][,1]) > 0),1]<-as.numeric(dBc[[i]][[j]][which(as.numeric(dBc[[i]][[j]][,1]) > 0),1])/maxP
        dBc[[i]][[j]][which(as.numeric(dBc[[i]][[j]][,1]) < 0),1]<-as.numeric(dBc[[i]][[j]][which(as.numeric(dBc[[i]][[j]][,1]) < 0),1])/maxN
        dBc[[i]][[j]][which(as.numeric(dBc[[i]][[j]][,1]) > 1),1]<-1
        dBc[[i]][[j]][which(as.numeric(dBc[[i]][[j]][,1]) < -1),1]<--1
      }
    }
  }
  dataOn<-new("GMMbyCond")
  dataBC(dataOn)<-dBc
  nO<-dBc
  for(l in 1:length(nO)){
    for(m in 1:length(nO[[l]])){
      dataIDon<-rownames(nO[[l]][[m]])[intersect(which(nO[[l]][[m]][,"clus"] == "P"), which(nO[[l]][[m]][,"status"] == "OK"))]
      sitesOn<-as.character(IDmap(dataGMM)[match(dataIDon, as.character(IDmap(dataGMM)[,"dataID"])),"S.cc"])
      nO[[l]][[m]]<-unique(as.character(sitesOn))
    }
  }
  speciesP(dataOn)<-nO
  return(dataOn)
}
