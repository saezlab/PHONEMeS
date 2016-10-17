#this is the script I used to get the data objects needed for PHONEMeS,
#based on the original objects that I was using before
#"~/Desktop/PHONEMeS_all/PHONEMeS_example/GMM_noFC.RData"
#"~/Desktop/PHONEMeS_all/PHONEMeS_example/dataObjects_2_2.RData"

# Working directory
# Set working directory to directory of this script (in RStudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# If you are using this file as Source, use:
# setwd(getSrcDirectory()[1])

# COMMENTS BY JAKOB
# Where is this file?
load("~/Desktop/PHONEMeS_all/PHONEMeS_example/GMM_noFC.RData")
colnames(GMM.res.ID)[6]<-"S.cc"
GMM.res.noFC<-GMM.res


load("../data/dataObjects_2_2.RData")

# Cannot be executed since GMM.res.noFC is missing
save(file="../data/dataObjects_PHONEMeS.RData", 
     list=c("GMM.res.ID","GMM.res.noFC", "GMM.res"))
