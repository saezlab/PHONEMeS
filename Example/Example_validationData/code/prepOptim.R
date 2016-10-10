#load the packages
library(BioNet)
library(igraph)
library(PHONEMeS)

# Set working directory to directory of this script (in RStudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# If you are using this file as Source, use:
# setwd(getSrcDirectory()[1])

#Load the network data and GMM results
load("../data/allD_noCSK_filt.RData")
load("../data/dataObjects_PHONEMeS.RData")
#Make the data objects that will be needed
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM.res, IDmap=GMM.res.ID, resFC=GMM.res.wFC)
#Choose the drug targets
targets.P<-list(cond1=c("MTOR_HUMAN"))
#Choose the drug treatments matching to the drug targets
#and match to what is present in the background network
data.P<-dataBycond(dataGMM, bg,scaled=TRUE,rowBycond=list(cond1=c("MTOR1 - Control2", "MTOR2 - Control2")))
show(data.P)
# Object of class GMMbyCond 
# GMM results list, by condition: 568 peptides 
# Conditions: 1 
# Species perturbed: 
#   cond1.MTOR1 - Control2 : 14 
#   cond1.MTOR2 - Control2 : 2
speciesP(data.P)
#Create the PKN list that will be used for optimisation
pknList<-buildNw(data.On=data.P, targets.On=targets.P, bg=bg,nK="no")
show(pknList)
# Object of class PKNlist 
# Interactions: 173 
# Direct interactions with targets: 4 
# Species: 107 
# Sinks: 13 
# Integrators: 31 
# Intermediates: 63 
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=3, nG=50)
save(file=paste("../cluster_scripts/data4cluster_23",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=4, nG=50)
save(file=paste("../cluster_scripts/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=5, nG=50)
save(file=paste("../cluster_scripts/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=6, nG=50)
save(file=paste("../cluster_scripts/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=7, nG=50)
save(file=paste("../cluster_scripts/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=8, nG=50)
save(file=paste("../cluster_scripts/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
