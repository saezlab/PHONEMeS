#load the packages
library(BioNet)
library(igraph)
library(PHONEMeS)
#Load the network data and GMM results
load("allD_noCSK_filt.RData")
load("dataObjects_PHONEMeS.RData")
#Make the data objects that will be needed
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM.res.noFC, IDmap=GMM.res.ID, resFC=GMM.res)
#Choose the drug targets
targets.P<-list(cond1=c("MTOR_HUMAN"))
#Choose the drug treatments matching to the drug targets
#and match to what is present in the background network
data.P<-dataBycond(dataGMM, bg,scaled=TRUE,rowBycond=list(cond1=c("MTOR1 - Control", "MTOR2 - Control")))
show(data.P)
speciesP(data.P)
#Create the PKN list that will be used for optimisation
pknList<-buildNw(data.On=data.P, targets.On=targets.P, bg=bg,nK="no")
show(pknList)
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=3, nG=50)
save(file=paste("data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=4, nG=50)
save(file=paste("data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=5, nG=50)
save(file=paste("data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=6, nG=50)
save(file=paste("data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=7, nG=50)
save(file=paste("data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=8, nG=50)
save(file=paste("ata4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
