buildNw<-function(data.On,targets.On, bg,
                    nK=c("all","no", "drugs2data", "data")){
      pkn<-buildPKN(data.On, targets.On, bg,nK=nK) 
      pknList<-PKNlist(pkn, targets.On, data.On)           
      return(pknList)   
}