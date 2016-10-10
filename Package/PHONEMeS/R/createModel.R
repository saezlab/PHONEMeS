createModel<-function(pknList=pknList,
                      AndBinF=rep(1, length(integrators(pknList))),
                      NoneBinF=rep(1, length(integrators(pknList))),
                      optParam=optParam){
  #complete.I.start is the same as above, for all
  complete.I<-interactions(pknList)
  sinks<-sinks(pknList)
  nodes<-species(pknList)
  complete.I.start<-interactionsD(pknList)
  intermediates<-intermediates(pknList)
  integrators<-integrators(pknList)
  #######create complete.I.sinks
  #sinks is always the same
  #nodes is always the same
  complete.I.sinks<-data.frame(S.AC=rep(NA, length(sinks)),
                               S.ID=rep(NA, length(sinks)),
                               K.AC=rep(NA, length(sinks)),
                               K.ID=rep(NA, length(sinks)),
                               res=rep(NA, length(sinks)),
                               pos=rep(NA, length(sinks)),
                               SID=rep(NA, length(sinks)),
                               S.cc=rep(NA, length(sinks)),
                               ntag=rep(NA, length(sinks)))
  sampling<-runif(n=length(sinks), 0,1)
  for(i in 1:length(sinks)){
    int<-complete.I[which(complete.I[,"S.cc"] == sinks[i]),]
    #now I need to copy through the interactions in int that have ntag>1
    ntagR<-do.call(rep,list(1:dim(int)[1], int$ntag))
    int<-int[ntagR,]
    n.nK<-length(grep(int[,"SID"], pattern="n"))
    n.nnK<-dim(int)[1] - n.nK 
    #this is the size of the bins
    #each non nK interac is worth 2 bins
    bins<-1/(n.nK + (2*n.nnK))
    boundary<-rep(2*bins, dim(int)[1])
    boundary[grep(int[,"SID"], pattern="n")]<-bins
    boundary<-cumsum(boundary)
    #this is the interaction that is chosen
    complete.I.sinks[i,]<-int[which(sampling[i] < boundary)[1],]
  }
  #######create complete.I.integrators
  complete.I.integrators<-complete.I[grep(complete.I[,"SID"],pattern="i"),]
  complete.I.integrators.tag<-rep(FALSE, dim(complete.I.integrators)[1])
  if(intgAsintm(optParam) == FALSE){
    #1 bin of the interval is for the AND (it gets corrected later), that's because 
    #it then gets less likely to be an AND as the number of inputs increases
    integrators.sampling<-runif(n=length(integrators), 0, 1)
    for(i in 1:length(integrators)){
      int<-complete.I.integrators[which(complete.I.integrators[,"S.cc"] == integrators[i]),]
      #now I need to copy through the interactions in int that have ntag>1
      ntagR<-do.call(rep,list(1:dim(int)[1], int$ntag))
      int<-int[ntagR,]
      #here I add 1 for the "nothing gets picked" bin
      bins<-1/(dim(int)[1]+AndBinF[i]+NoneBinF[i])
      boundary<-rep(bins, dim(int)[1])
      boundary<-cumsum(boundary)
      #if we are above, then it's an AND
      if(all(boundary < integrators.sampling[i])){
        if(integrators.sampling[i] <  ((dim(int)[1]+AndBinF[i]) * bins)){
          complete.I.integrators.tag[match(int$SID, complete.I.integrators$SID)]<-TRUE
        }
        #if this is not the case then no integrator is picked and nothing gets turned to TRUE
      }else{
        #we pick the interaction depending on where we fall
        complete.I.integrators.tag[match(int$SID[which(integrators.sampling[i] < boundary)[1]], complete.I.integrators$SID)]<-TRUE
        #this bit will be used when we have weights to correct, in the first step
        #when there is only one in edge it always gets in
        #if we're below and there was only 1 interaction, then it's actually a FALSE
        #this allows me to then treat the one incoming edge as and AND
      }
    }
  }else{
    for(i in 1:length(integrators)){
      int<-complete.I.integrators[which(complete.I.integrators[,"S.cc"] == integrators[i]),]
      #now I need to copy through the interactions in int that have ntag>1
      ntagR<-do.call(rep,list(1:dim(int)[1], int$ntag))
      int<-int[ntagR,]
      sampling<-runif(n=dim(int)[1], 0,1)
      bins<-1/dim(int)[1]
      sI<-which(sampling < bins)
      int.sI<-unique(int[sI,"SID"])
      complete.I.integrators.tag[which(complete.I.integrators[,"SID"] %in% int.sI)]<-TRUE
      }
  }
  complete.I.integrators<-complete.I.integrators[complete.I.integrators.tag,]
  ########create complete.I.intermediates
  complete.I.intermediates<-complete.I[(complete.I[,"S.cc"] %in% intermediates),]
  #intermediates.t is the same for all
  complete.I.intermediates.tag<-rep(FALSE, dim(complete.I.intermediates)[1])
  for(i in 1:length(intermediates)){
    int<-complete.I.intermediates[which(complete.I.intermediates[,"S.cc"] == intermediates[i]),]
    #now I need to copy through the interactions in int that have ntag>1
    ntagR<-do.call(rep,list(1:dim(int)[1], int$ntag))
    int<-int[ntagR,]
    sampling<-runif(n=dim(int)[1], 0,1)
    n.nK<-length(grep(int[,"SID"], pattern="n"))
    n.nnK<-dim(int)[1] - n.nK
    bins<-1/(n.nK + (2*n.nnK))
    ps<-rep(2*bins, dim(int)[1])
    ps[grep(int[,"SID"], pattern="n")]<-bins
    sI<-which(sampling < ps)
    int.sI<-unique(int[sI,"SID"])
    complete.I.intermediates.tag[which(complete.I.intermediates[,"SID"] %in% int.sI)]<-TRUE
  }
  complete.I.intermediates<-complete.I.intermediates[complete.I.intermediates.tag,]
  if(cstart(optParam) == TRUE){
    model<-rbind(complete.I.start, complete.I.integrators, complete.I.sinks, complete.I.intermediates)
  }else{
    model<-rbind(complete.I.integrators, complete.I.sinks, complete.I.intermediates)
  }
  return(model)
}