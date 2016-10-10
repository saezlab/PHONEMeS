pathsD<-function(nwTable, nM,nodes2link,drug){
#this makes the weighted network object that we need
ge <- new("graphNEL", nodes = unique(c(nwTable[,"K.ID"], nwTable[,"S.cc"])), edgemode = "directed")
ep<-nwTable$ntag
ep[which(is.na(ep))]<-0
ep<-ep[match(names(ep), nwTable$SID)]
complete.w <- addEdge(nwTable[,"K.ID"], nwTable[,"S.cc"], ge, weights=(nM-ep))
complete.w<-igraph.from.graphNEL(complete.w)
#this is the nodes that we want to connect
dataVn<-match(nodes2link, V(complete.w)$name)
dataVn<-dataVn[!is.na(dataVn)]
#this is the kinase that we want to connect the nodes to
targetVn<-match(drug, V(complete.w)$name)
Paths<-get.all.shortest.paths(g=complete.w, from=targetVn, to=dataVn, mode="out", weights=NULL)
#this is the number of nodes for which a path was found
#length(Paths$res)
#this is the original number of nodes
#length(dataVn)
print(paste("Paths from your drug target to your data nodes were found for ", length(Paths$res), " nodes", sep=""))
#this maps the paths back to the original edge table
tp1<-V(complete.w)$name[Paths$res[[1]]]
DrugPaths<-nwTable[intersect(which(nwTable$K.ID == tp1[1]),which(nwTable$S.cc == tp1[2])),]
for(j in 1:(length(tp1)-1)){
    temp<-nwTable[intersect(which(nwTable$K.ID == tp1[j]),which(nwTable$S.cc == tp1[j+1])),]
    temp<-temp[which(temp$ntag == max(temp$ntag)),]
    DrugPaths<-rbind(DrugPaths, temp[1,])
  }
for(i in 2:length(Paths$res)){
  tp1<-V(complete.w)$name[Paths$res[[i]]]
  for(j in 1:(length(tp1)-1)){
    temp<-nwTable[intersect(which(nwTable$K.ID == tp1[j]),which(nwTable$S.cc == tp1[j+1])),]
    temp<-temp[which(temp$ntag == max(temp$ntag)),]
    DrugPaths<-rbind(DrugPaths, temp[1,])
  }
}
return(DrugPaths)
}
