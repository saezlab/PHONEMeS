integratorWeights<-function(models.sel=models.sel, 
                            pknList=pknList,
                            models.sel.intgAnd=models.sel.intgAnd,
                            models.sel.comb=models.sel.comb,
                            optParam=optParam){
  integrators<-integrators(pknList)
  complete.I<-interactions(pknList)                       
  integrators2paste<-data.frame(S.AC=rep(NA, length(integrators(pknList))*cap(optParam)),
                                S.ID=rep(NA, length(integrators(pknList))*cap(optParam)),
                                K.AC=rep(NA, length(integrators(pknList))*cap(optParam)),
                                K.ID=rep(NA, length(integrators(pknList))*cap(optParam)),
                                res=rep(NA, length(integrators(pknList))*cap(optParam)),
                                pos=rep(NA, length(integrators(pknList))*cap(optParam)),
                                SID=rep(NA, length(integrators(pknList))*cap(optParam)),
                                S.cc=rep(NA, length(integrators(pknList))*cap(optParam)),
                                ntag=rep(NA, length(integrators(pknList))*cap(optParam)))
  #AND bin factor if intgAsintm=FALSE:it is a single number that adds to the +1 in the
  #determination of the size of bins 
  #I did a simple proportional, same as for the others, that gets added up
  #this effectively increases the number of AND bins, by reducing
  #the size of all bins and allocating more to AND
  if(intgAsintm(optParam) == FALSE){  
    whichI<-rep(NA, length(models.sel))
    AndBinF<-rep(1, length(integrators(pknList)))
    NoneBinF<-rep(1, length(integrators(pknList)))
    count<-1
    for(i in 1:length(integrators(pknList))){
      AndCount<-0
      NoneCount<-0
      for(m in 1:length(models.sel)){
        tempI<-models.sel[[m]][which(models.sel[[m]][,"S.cc"] == integrators(pknList)[i]),]
        if(dim(tempI)[1] == 1){
          whichI[m]<-tempI[,"SID"]
        }else{
          if(dim(tempI)[1] > 1){
            AndCount<-AndCount+1
          }else{
            NoneCount<-NoneCount+1
          }
        }
      }
      NoneBinF[i]<-cap(optParam)*(NoneCount/length(models.sel))
      NoneBinF[i]<-ceiling(NoneBinF[i])	
      if(!all(is.na(whichI))){
        AndBinF[i]<-cap(optParam)*(AndCount/length(models.sel))
        AndBinF[i]<-ceiling(AndBinF[i])
        #now I need a table of stuff to copy - for this I can do the same as for sinks
        occ<-table(whichI)
        na<-sum(is.na(whichI))
        if(na > 0){
          occ<-c(occ, na)
          names(occ)[length(occ)]<-"na"
        }
        freq<-round((occ/sum(occ))*cap(optParam), digits=0)
        if(length(freq) == 1){
          integrators2paste[count:(count+cap(optParam)-1),]<-interactions(pknList)[match(names(freq)[1], interactions(pknList)[,"SID"]),]
          count<-count+cap(optParam)
        }else{
          if(na > 0){
            freq<-freq[-which(names(freq) == "na")]
          }
          for(j in 1:length(freq)){
            if(freq[j] != 0){
              integrators2paste[count:(count+freq[j]-1),]<-interactions(pknList)[match(names(freq)[j], interactions(pknList)[,"SID"]),]
              count<-count+freq[j]
            }   	
          }
        }
      }
    }
    #this is the bit where intg are treated as intm, if the param is true
  }else{
    AndBinF<-colSums(models.sel.intgAnd)/length(models.sel)
    #this becomes the new flip boundary
    #then for the copy of the interactions I can do the same as sinks
    count<-1
    for(s in 1:length(integrators(pknList))){
      tempintM<-models.sel.comb[which(models.sel.comb[,"S.cc"] == integrators(pknList)[s]),]
      occ<-table(tempintM[,"SID"])
      na<-sum(is.na(tempintM[,"SID"]))
      if(na > 0){
        occ<-c(occ, na)
        names(occ)[length(occ)]<-"na"
      }
      freq<-round((occ/sum(occ))*cap(optParam), digits=0)
      if(length(freq) != 0){
        if(length(freq) == 1){
          integrators2paste[count:(count+cap(optParam)-1),]<-tempintM[1,]
          count<-count+cap(optParam)
        }else{
          if(na > 0){
            freq<-freq[-which(names(freq) == "na")]
          }
          for(j in 1:length(freq)){
            if(freq[j] != 0){
              integrators2paste[count:(count+freq[j]-1),]<-tempintM[match(names(freq)[j], tempintM[,"SID"]),]
              count<-count+freq[j]
            }   
          }
        }
      }    
    }	
  }
  integrators2paste<-integrators2paste[!is.na(integrators2paste[,"S.cc"]),]
  return(list(integrators2paste=integrators2paste, AndBinF=AndBinF, NoneBinF=NoneBinF))
}