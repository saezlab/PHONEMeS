#this is the script I used to get the data objects needed for PHONEMeS,
#based on the original objects that I was using before
#"~/Desktop/PHONEMeS_all/PHONEMeS_example/GMM_noFC.RData"
#"~/Desktop/PHONEMeS_all/PHONEMeS_example/dataObjects_2_2.RData"
load("~/Desktop/PHONEMeS_all/PHONEMeS_example/GMM_noFC.RData")
colnames(GMM.res.ID)[6]<-"S.cc"
GMM.res.noFC<-GMM.res
load("~/Desktop/PHONEMeS_all/PHONEMeS_example/dataObjects_2_2.RData")
save(file="~/Desktop/PHONEMeS_all/S4/dataObjects_PHONEMeS.RData", 
     list=c("GMM.res.ID","GMM.res.noFC", "GMM.res"))
