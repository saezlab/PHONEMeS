buildPKN <-function(data.On,targets.On, bg,
                    nK=c("all","no", "drugs2data", "data")){
  dataNodes<-unique(unlist(speciesP(data.On)))
  drugNodes<-unique(unlist(targets.On))
  pSdf<-interactions(bg)
  #this takes care of the case where we don't want any nK interacs
  if(nK == "no"){
    pSdf<-pSdf[-grep("n", pSdf$SID),]
  }		
  #make a graph object from the network - at the protein level
  allD.na<-union(which(is.na(pSdf[,"K.ID"])), which(is.na(pSdf[,"S.ID"])))
  allD.nw<-makeNetwork(source=pSdf[-allD.na,"K.ID"], 
                       target=pSdf[-allD.na,"S.ID"], 
                       edgemode="directed")
  allD.nw<-igraph.from.graphNEL(allD.nw)
  #these are needed if we want nK interacs in some parts of the network
  if(nK == "drugs2data" || nK == "data"){
    #make a version of pSdf that has no NK interacs
    pSdf.noN<-pSdf[-grep("n", pSdf$SID),]
    allD.na<-pSdf[-allD.na,]
    allD.na.noN<-allD.na[-grep("n", allD.na$SID),]
    allD.nw.noN<-makeNetwork(source=allD.na.noN[,"K.ID"], 
                             target=allD.na.noN[,"S.ID"], 
                             edgemode="directed")
    allD.nw.noN<-igraph.from.graphNEL(allD.nw.noN)
  }	
  #########Extract the interactions linking the drug targets to each other
  if(nK == "drugs2data" || nK == "data"){
    nw.2use<-allD.nw.noN
    pSdf.2use<-pSdf.noN
  }else{
    nw.2use<-allD.nw
    pSdf.2use<-pSdf
  }
  #
  if(length(drugNodes) != 1){
    targetVn<-match(drugNodes, V(nw.2use)$name)
    #for each drug target, look for all shortest paths to all other targets (max 7 nodes) 
    targetPaths<-rep(NA,2)
    for(t in 1:length(targetVn)){
      Paths<-get.all.shortest.paths(g=nw.2use, 
                                    from=targetVn[t], 
                                    to=targetVn[-t], 
                                    mode="out", weights=NA)
      Paths<-lapply(Paths$res, function(x){return(V(nw.2use)$name[x])})
      Paths.l<-lapply(Paths, length)
      Paths.l<-unlist(Paths.l)
      Paths<-Paths[which(Paths.l <= 7)]
      Paths<-lapply(Paths, function(x){if(length(x) == 2){return(x)};
                                       if(length(x) == 3){return(rbind(x[1:2], x[2:3]))};
                                       if(length(x) == 4){
                                         return(rbind(x[1:2], x[2:3], x[3:4]))
                                       };
                                       if(length(x) == 6){
                                         return(rbind(x[1:2], x[2:3], x[3:4], x[4:5], x[5:6]))
                                       };
                                       if(length(x) == 7){
                                         return(rbind(x[1:2], x[2:3], x[3:4], x[4:5], x[5:6], x[6:7]))
                                       };
                                       if(length(x) == 5){
                                         return(rbind(x[1:2], x[2:3], x[3:4], x[4:5]))}
      })
      if(length(Paths) != 0){
        if(length(Paths) > 1){
          for(i in 2:length(Paths)){
            Paths[[i]]<-rbind(Paths[[i-1]], Paths[[i]])
          }
        }
        targetPaths<-rbind(targetPaths, Paths[[length(Paths)]]) 
      }
    }
    if(!all(is.na(targetPaths))){
      #remove the initialiser line and the duplicates
      targetPaths<-targetPaths[2:dim(targetPaths)[1],]
      targetPaths<-targetPaths[!duplicated(targetPaths),]
      print(paste("The following drug targets could not be included in the drug targets network:", 
                  drugNodes[!(drugNodes %in% c(targetPaths[,1], targetPaths[,2]))], sep=" "))
      #having obtained the protein level paths, I now look for all underlying interactions
      allD.2keep<-apply(pSdf.2use, MARGIN=1, 
                        function(x){k<-which(targetPaths[,1] == x["K.ID"]);
                                    s<-which(targetPaths[,2] == x["S.ID"]);
                                    return(length(intersect(k,s) != 0))})
      targetsPath.I<-pSdf.2use[which(allD.2keep == 1),]	
      print(paste("The network connecting your drug targets contains", 
                  dim(targetsPath.I)[1], "interactions and", 
                  length(unique(c(targetsPath.I[,"S.ID"], targetsPath.I[,"K.ID"]))),"proteins", sep=" "))
      #this is the bit that adds the drugT that cannot be connected with each other
      #but that still do appear in the network
      print(paste("The following drug targets could not be included:", 
                  drugNodes[!(drugNodes %in% pSdf.2use[,"K.ID"])], sep=" "))	
      dN.off<-drugNodes[!(drugNodes %in% c(targetPaths[,1], targetPaths[,2]))]
      dN.in<-dN.off[(dN.off %in% pSdf.2use[,"K.ID"])]
      if(length(dN.in) != 0){
        temp<-pSdf[1:length(dN.in),]
        temp[,1:dim(temp)[2]]<-NA
        temp[,c("S.ID","K.ID")]<-cbind(dN.in, dN.in)
        targetsPath.I<-rbind(targetsPath.I, temp)
      }
    }else{
      #this is the bit that happens if the drug targets cannot be connected
      print(paste("The following drug targets could not be included:", 
                  drugNodes[!(drugNodes %in% pSdf.2use[,"K.ID"])], sep=" "))
      dN.in<-drugNodes[(drugNodes %in% pSdf.2use[,"K.ID"])]
      if(length(dN.in) != 0){
        targetsPath.I<-pSdf[1:length(dN.in),]
        targetsPath.I[,1:dim(targetsPath.I)[2]]<-NA
        targetsPath.I[,c("S.ID","K.ID")]<-cbind(dN.in, dN.in)
      }
    }	
    #this is what happens if there is only 1 drug target
  }else{
    targetsPath.I<-pSdf[1,]
    targetsPath.I[1,1:dim(targetsPath.I)[2]]<-NA
    targetsPath.I[,c("S.ID","K.ID")]<-drugNodes
  }
  ######Extract the interactions linking data sites to their upstream kinase
  allD.2keep<-pSdf[,"S.cc"] %in% dataNodes
  dataSites.I<-pSdf[which(allD.2keep == TRUE),]
  print(paste("Your data contains information about", length(unique(dataNodes)), 
              "sites, of which", length(unique(dataSites.I[,"S.cc"])), 
              "are in your network", sep=" "))
  ######Extract the paths from data network to targets network: 
  if(nK == "data"){
    nw.2use<-allD.nw.noN
    pSdf.2use<-pSdf.noN
  }else{
    nw.2use<-allD.nw
    pSdf.2use<-pSdf
  }
  #link targetsPath.I[,"S.ID"] (or targetsPath.I[,"K.ID"]) to dataSites.I[,"K.ID"]
  dataKVn<-match(unique(dataSites.I[,"K.ID"]), V(nw.2use)$name)
  targetsSVn<-match(unique(c(targetsPath.I[,"S.ID"], targetsPath.I[,"K.ID"])), 
                    V(nw.2use)$name)
  #I don't want to match the stuff in dataKVn that is already in targetsSVn
  dataKVn<-dataKVn[!(dataKVn %in% targetsSVn)]
  dataKVn<-dataKVn[!is.na(dataKVn)]
  data2targetPaths<-rep(NA,2)
  for(t in 1:length(dataKVn)){
    Paths<-get.all.shortest.paths(g=nw.2use, 
                                  from=dataKVn[t], to=targetsSVn, 
                                  mode="in", weights=NA)
    Paths<-lapply(Paths$res, function(x){return(V(nw.2use)$name[x])})
    Paths<-lapply(Paths, function(x){if(length(x) == 2){return(x)};
                                     if(length(x) == 3){return(rbind(x[2:1], x[3:2]))};
                                     if(length(x) == 4){return(rbind(x[2:1], x[3:2], x[4:3]))};
                                     if(length(x) >= 5){return(NA)}})
    if(length(Paths) != 0){
      if(length(Paths) == 1){
        data2targetPaths<-rbind(data2targetPaths, Paths[[length(Paths)]])
      }else{
        for(i in 2:length(Paths)){
          if(length(Paths[[i]]) == 1){
            Paths[[i]]<-Paths[[i-1]]
          }else{
            Paths[[i]]<-rbind(Paths[[i-1]], Paths[[i]])
          }
        }
        data2targetPaths<-rbind(data2targetPaths, Paths[[length(Paths)]]) 
      }  
    }
  }
  #remove the initialiser line and the duplicates
  data2targetPaths<-data2targetPaths[2:dim(data2targetPaths)[1],]
  data2targetPaths<-data2targetPaths[!duplicated(data2targetPaths),]
  #these paths are at the protein level: add underlying K-S interactions
  allD.2keep<-apply(data2targetPaths, MARGIN=1, 
                    function(x){k<-which(pSdf.2use[,"K.ID"] == x[1]);
                                s<-which(pSdf.2use[,"S.ID"] == x[2]);
                                return(intersect(k,s))})
  data2targetPaths.I<-pSdf.2use[unique(unlist(allD.2keep)),]
  if(length(drugNodes) != 1){
    complete.I<-rbind(dataSites.I, targetsPath.I, data2targetPaths.I)
  }else{
    complete.I<-rbind(dataSites.I, data2targetPaths.I)
  }
  complete.I<-complete.I[!duplicated(complete.I),]
  print(paste("Your complete network contains",dim(complete.I)[1],
              "kinase/phosphatase substrate interactions (possibly non unique)", sep=" "))
  ######Add the integrator information
  #for each protein that appears both in substrate and in kinase, I need to create a
  #link from the substrate node (with its site), to the kinase node
  toadd<-data.frame(
    complete.I[which(complete.I[,"S.ID"] %in% complete.I[,"K.ID"]),"S.cc"],
    complete.I[which(complete.I[,"S.ID"] %in% complete.I[,"K.ID"]),"S.ID"],
    stringsAsFactors=FALSE)
  toadd<-toadd[!duplicated(toadd),]
  colnames(toadd)<-c("K.ID", "S.cc")
  toadd$S.AC<-rep(NA, dim(toadd)[1])
  toadd$S.ID<-rep(NA, dim(toadd)[1])
  toadd$K.AC<-rep(NA, dim(toadd)[1])
  toadd$res<-rep(NA, dim(toadd)[1])
  toadd$pos<-rep(NA, dim(toadd)[1])
  toadd$SID<-paste("i", seq(1,dim(toadd)[1]), sep="")
  complete.I<-rbind(complete.I, toadd)
  print(paste(dim(toadd)[1], "interactions to",length(unique(toadd[,"S.cc"])),
              "integrator nodes were added", sep=" "))
  complete.I$ntag<-rep(1, dim(complete.I)[1])
  if(length(which(is.na(complete.I$K.ID))) != 0){
    complete.I<-complete.I[-which(is.na(complete.I$K.ID)),]
  }
  if(length(which(is.na(complete.I$SID))) != 0){
    complete.I<-complete.I[-which(is.na(complete.I$SID)),]
  }
  PKN<-new("KPSbg", interactions=complete.I, species=unique(c(complete.I$K.ID, complete.I$S.cc)))
  return(PKN)
}
