# Set working directory to directory of this script (in RStudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# If you are using this file as Source, use:
# setwd(getSrcDirectory()[1])

## COMMENT BY JAKOB:
# Where is this file? - This script does not really make sense without the original .RData file...
# Anyway, I will adjust the output paths for consistency
load(file="~/Desktop/Ed-data/data2.RData")
#
sum(duplicated(data2$psite))#906
sum(duplicated(data2))#[1] 0
sum(duplicated(data2$db.ID))#0
#there are duplicated peptides...
#they are not duplicated row though so I assume there is a difference
dup1<-data2[duplicated(data2$psite),]
data2[which(data2$psite == dup1$psite[1]),]
#there seems to be one with a better delta score and lower pFDR, that also
#has less missing values - maybe I should pick the one with the lowest pFDR?
row2keep<-!duplicated(data2$psite)
n0s<-apply(data2[,32:175], MARGIN=1, function(x){sum(is.na(x))})
for(i in 1:dim(data2)[1]){
  if(row2keep[i] == FALSE){
    dr<-data2[which(data2$psite == data2$psite[i]),"pFDR"]
    sc<-data2[which(data2$psite == data2$psite[i]),"max_delta_score"]
    nm<-n0s[which(data2$psite == data2$psite[i])]
    if(data2[i,"pFDR"] == min(dr) && data2[i,"max_delta_score"] == max(sc) && n0s[i] == min(nm)) row2keep[i]<-TRUE
  }
}
sum(row2keep)#[1] 12694
sum(duplicated(data2$psite[row2keep]))#195
#there are still some left
#that's because of the NA min condition
#the stuff below does the same as above but removing the NA condition
data2.noDup<-data2
data2.noDup$psite<-rep(NA, length(data2.noDup$psite))
for(i in 1:dim(data2)[1]){
  dr<-data2[which(data2$psite == data2$psite[i]),]
  if(dim(dr)[1] == 1){
    data2.noDup[i,"psite"]<-data2[i,"psite"]
  }else{
    if((data2[i,"pFDR"] == min(dr$pFDR)) && (data2[i,"max_delta_score"] == max(dr$max_delta_score))) {
      data2.noDup[i,"psite"]<-data2[i,"psite"]
    }
  }
}
#
sum(is.na(data2.noDup$psite)) #1139
sum(duplicated(data2$psite[!is.na(data2.noDup$psite)]))
d12<-data2[!is.na(data2.noDup$psite),32:175]
rownames(d12)<-data2[!is.na(data2.noDup$psite),"psite"]
save(d12, file="../data/MSdata.RData")
