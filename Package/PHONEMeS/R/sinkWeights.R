sinkWeights <-
  function(pknList=pknList,
           optParam=optParam,
           models.sel.comb=models.sel.comb, 
           models.sel=models.sel){
    sinks<-sinks(pknList)
    cap<-cap(optParam)
    sink2paste<-data.frame(S.AC=rep(NA, length(sinks)*cap),
                           S.ID=rep(NA, length(sinks)*cap),
                           K.AC=rep(NA, length(sinks)*cap),
                           K.ID=rep(NA, length(sinks)*cap),
                           res=rep(NA, length(sinks)*cap),
                           pos=rep(NA, length(sinks)*cap),
                           SID=rep(NA, length(sinks)*cap),
                           S.cc=rep(NA, length(sinks)*cap),
                           ntag=rep(NA, length(sinks)*cap))
    count<-1
    for(s in 1:length(sinks)){
      tempsinkM<-models.sel.comb[which(models.sel.comb[,"S.cc"] == sinks[s]),]
      occ<-table(tempsinkM[,"SID"])
      if(any(occ == length(models.sel)) && length(occ) != 1){
        occ[which(occ == length(models.sel))]<-occ[order(occ, decreasing=TRUE)[2]]
        freq<-round((occ/sum(occ))*cap, digits=0)
      }else{
        freq<-round((occ/sum(occ))*cap, digits=0)
      }
      if(length(freq) == 1){
        sink2paste[count:(count+cap-1),]<-tempsinkM[1,]
        count<-count+cap
      }else{
        for(j in 1:length(freq)){
          if(freq[j] != 0){
            sink2paste[count:(count+freq[j]-1),]<-tempsinkM[match(names(freq)[j], tempsinkM[,"SID"]),]
            count<-count+freq[j]
          }   
        }
      }
    }
    sink2paste<-sink2paste[!is.na(sink2paste[,"S.cc"]),]
    return(sink2paste)
  }
