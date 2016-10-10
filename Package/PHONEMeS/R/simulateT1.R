simulateT1 <-
  function(nodes.t=nodes.t,
           intermediates.AND=intermediates.AND,
           integrators.AND=integrators.AND,
           intm.inN=intm.inN,
           intg.inN=intg.inN,
           sink.inN=sink.inN){
    for(i in 1:length(intm.inN)){
      nodes.t[names(intm.inN)[i],2]<-ifelse(intermediates.AND[i],min(nodes.t[intm.inN[[i]],1]),max(nodes.t[intm.inN[[i]],1]))
    }
    if(all(is.na(integrators.AND))){
      for(i in 1:length(intg.inN)){
        nodes.t[names(intg.inN)[i],2]<-min(nodes.t[intg.inN[[i]],1])
      }
    }else{
      for(i in 1:length(intg.inN)){
        nodes.t[names(intg.inN)[i],2]<-ifelse(integrators.AND[i],min(nodes.t[intg.inN[[i]],1]),max(nodes.t[intg.inN[[i]],1]))
      }
    }	
    for(i in 1:length(sink.inN)){
      nodes.t[names(sink.inN)[i],2]<-max(nodes.t[sink.inN[[i]],1])
    }
    return(nodes.t)
  }
