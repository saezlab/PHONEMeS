intermediateWeights <-
  function(pknList=pknList,
           models.sel.intAnd=models.sel.intAnd,
           models.sel=models.sel,
           models.sel.comb=models.sel.comb,
           optParam=optParam){
    AndFlipF<-colSums(models.sel.intAnd)/length(models.sel)
    #this becomes the new flip boundary
    #then for the copy of the interactions I can do the same as sinks
    interm2paste<-data.frame(S.AC=rep(NA, length(intermediates(pknList))*cap(optParam)),
                             S.ID=rep(NA, length(intermediates(pknList))*cap(optParam)),
                             K.AC=rep(NA, length(intermediates(pknList))*cap(optParam)),
                             K.ID=rep(NA, length(intermediates(pknList))*cap(optParam)),
                             res=rep(NA, length(intermediates(pknList))*cap(optParam)),
                             pos=rep(NA, length(intermediates(pknList))*cap(optParam)),
                             SID=rep(NA, length(intermediates(pknList))*cap(optParam)),
                             S.cc=rep(NA, length(intermediates(pknList))*cap(optParam)),
                             ntag=rep(NA, length(intermediates(pknList))*cap(optParam)))
    count<-1
    for(s in 1:length(intermediates(pknList))){
      tempintM<-models.sel.comb[which(models.sel.comb[,"S.cc"] == intermediates(pknList)[s]),]
      occ<-table(tempintM[,"SID"])
      na<-sum(is.na(tempintM[,"SID"]))
      if(na > 0){
        occ<-c(occ, na)
        names(occ)[length(occ)]<-"na"
      }
      freq<-round((occ/sum(occ))*cap(optParam), digits=0)
      if(length(freq) != 0){
        if(length(freq) == 1){
          interm2paste[count:(count+cap(optParam)-1),]<-tempintM[1,]
          count<-count+cap(optParam)
        }else{
          for(j in 1:length(freq)){
            if(na > 0){
              freq<-freq[-which(names(freq) == "na")]
            }
            if(freq[j] != 0){
              interm2paste[count:(count+freq[j]-1),]<-tempintM[match(names(freq)[j], tempintM[,"SID"]),]
              count<-count+freq[j]
            }   
          }
        }
      }    
    }
    interm2paste<-interm2paste[!is.na(interm2paste[,"S.cc"]),]
    return(list(interm2paste=interm2paste, AndFlipF=AndFlipF))
  }
