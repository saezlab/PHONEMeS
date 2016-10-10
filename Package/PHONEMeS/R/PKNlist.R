PKNlist<-function(PKN, targets.On, data.On){
  complete.I<-interactions(PKN)
  dataNodes<-unique(unlist(speciesP(data.On)))
  drugNodes<-unique(unlist(targets.On))
  #combine duplicate edges and adjust the weights accordingly
  complete.I.2<-complete.I
  ntag.2<-complete.I[,9]
  tag<-rep(TRUE, dim(complete.I.2)[1])
  for(i in 1:dim(complete.I.2)[1]){
    w<-intersect(which(complete.I.2[,"K.ID"] == complete.I.2[i,"K.ID"]),
                 which(complete.I.2[,"S.cc"] == complete.I.2[i,"S.cc"]))
    if(all(tag[w] == TRUE)){
      if(length(w) > 1){
        ntag.2[i]<-sum(ntag.2[w])
        tag[w]<-FALSE
        tag[i]<-TRUE
        complete.I.2[i,"SID"]<-paste(complete.I[w,"SID"], collapse=";")
      }
    }
  }
  complete.I.2[,"ntag"]<-ntag.2
  complete.I.2<-complete.I.2[tag,]
  #get the start network
  complete.I.start<-complete.I.2[intersect(which(complete.I.2[,"K.ID"] %in% drugNodes), 
                                           which(complete.I.2[,"S.cc"] %in% dataNodes)),]
  if(length(grep(complete.I.start[,"SID"], pattern="n")) != 0){
    whichN<-grep(complete.I.start[,"SID"], pattern="n")
    whichN<-whichN[!(whichN %in% c(grep(complete.I.start[,"SID"], pattern="e"), grep(complete.I.start[,"SID"], pattern="p"),
                                   grep(complete.I.start[,"SID"], pattern="h"),grep(complete.I.start[,"SID"], pattern="d")))]
    complete.I.start<-complete.I.start[-whichN,]
  }
  #there are still edges that reach a sink but then do not go up...
  #they should be removed
  ksink<-complete.I.2[which(!(complete.I.2$S.cc %in% complete.I.2$K.ID)),"K.ID"]
  ksink<-ksink[which(!ksink %in% complete.I.2$S.ID)]
  if(any(ksink %in% drugNodes)) ksink<-ksink[which(!ksink %in% drugNodes)]
  complete.I.2<-complete.I.2[-which(complete.I.2$K.ID %in% ksink), ]
  nodes<-unique(c(complete.I.2[,"K.ID"], complete.I.2[,"S.cc"]))
  nodes<-nodes[!is.na(nodes)]
  sinks<-nodes[!(nodes %in% complete.I.2[,"K.ID"])]
  integrators<-unique(complete.I.2[grep(complete.I.2[,"SID"],pattern="i"),"S.cc"])
  integrators<-integrators[!is.na(integrators)]
  intermediates<-setdiff(nodes, sinks)
  intermediates<-setdiff(intermediates, integrators)
  #this makes sure that we only look at intermediates that have an incoming edge
  intermediates<-intermediates[intermediates %in% complete.I.2[,"S.cc"]]
  PKN.list<-new("PKNlist", interactions=complete.I.2,
                interactionsD=complete.I.start,
                species=nodes,
                sinks=sinks,
                integrators=integrators,
                intermediates=intermediates)
  return(PKN.list)
}