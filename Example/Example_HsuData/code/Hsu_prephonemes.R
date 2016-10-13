
# Set working directory to directory of this script (in RStudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# If you are using this file as Source, use:
# setwd(getSrcDirectory()[1])

#load the packages
library(BioNet)
library(igraph)
library(PHONEMeS)
#Load the network data and GMM results
load(file="../data/data_Hsu_PHONEMeS_Feb13.RData")
load("../../Example_MainData/data/allD_noCSK_filt.RData")
#Make the data objects that will be needed
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM, IDmap=GMM.ID, resFC=GMM.wFC)
#Choose the drug targets
targets.P<-list(cond1=c("MTOR_HUMAN"))
#Choose the drug treatments matching to the drug targets
#and match to what is present in the background network
data.P<-dataBycond(dataGMM, bg,scaled=TRUE,rowBycond=list(cond1=c("rapa - ins")))
show(data.P)
#GMM results list, by condition: 603 peptides 
#Conditions: 1 
#Species perturbed: 
#  cond1.rapa - ins : 38 
speciesP(data.P)
#$cond1
#$cond1$`rapa - ins`
#[1] "PDCD4_HUMAN.S.76"   "4EBP2_HUMAN.Y.34"   "PDCD4_HUMAN.S.457"  "KIF4A_HUMAN.S.801"  "NPM_HUMAN.S.70"    
#[6] "KS6B1_HUMAN.S.447"  "TBCD4_HUMAN.T.752"  "RPTOR_HUMAN.S.863"  "PYR1_HUMAN.S.1859"  "AKT3_HUMAN.S.472"  
#[11] "NDRG1_HUMAN.S.332"  "TPIS_HUMAN.S.21"    "P53_HUMAN.S.315"    "PRP4B_HUMAN.S.578"  "HS90B_HUMAN.S.226" 
#[16] "NSUN2_HUMAN.S.743"  "SRRM1_HUMAN.S.738"  "SRRM2_HUMAN.S.2449" "VINEX_HUMAN.S.188"  "HDAC2_HUMAN.S.516" 
#[21] "XPO6_HUMAN.T.204"   "APC1_HUMAN.S.688"   "FNBP1_HUMAN.S.296"  "ADDG_HUMAN.S.683"   "ELAV1_HUMAN.S.202" 
#[26] "PRP4B_HUMAN.S.277"  "BICD2_HUMAN.S.582"  "BORG5_HUMAN.S.350"  "ATX2_HUMAN.T.666"   "SRRM2_HUMAN.S.1101"
#[31] "SCRIB_HUMAN.S.1486" "NU153_HUMAN.S.334"  "MYCN_HUMAN.T.58"    "TP53B_HUMAN.S.1462" "DNLI1_HUMAN.T.233" 
#[36] "F122A_HUMAN.S.143"  "SRRM1_HUMAN.T.581"  "TP53B_HUMAN.S.1426"
#Create the PKN list that will be used for optimisation
pknList<-buildNw(data.On=data.P, targets.On=targets.P, bg=bg,nK="no")
#[1] "Your data contains information about 38 sites, of which 15 are in your network"
#[1] "Your complete network contains 156 kinase/phosphatase substrate interactions (possibly non unique)"
#[1] "40 interactions to 20 integrator nodes were added"
show(pknList)
#Interactions: 109 
#Direct interactions with targets: 2 
#Species: 73 
#Sinks: 13 
#Integrators: 20 
#Intermediates: 40
#To have a look at what your data contains: show(data.P)
#To have a look at the sites perturbed under each drug: speciesP(data.P)
#To have a look at your network: show(pknList)
#set the parameters for optimisation
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=9, nG=50)
save(file=paste("../cluster_scripts/rapa/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=10, nG=50)
save(file=paste("../cluster_scripts/rapa/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=11, nG=50)
save(file=paste("../cluster_scripts/rapa/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=12, nG=50)
save(file=paste("../cluster_scripts/rapa/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=13, nG=50)
save(file=paste("../cluster_scripts/rapa/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=14, nG=50)
save(file=paste("../cluster_scripts/rapa/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
#############################################
#Torin
#load the packages
library(BioNet)
library(igraph)
library(PHONEMeS)
#Load the network data and GMM results
load(file="../data/data_Hsu_PHONEMeS_Feb13.RData")
load("../../Example_MainData/data/allD_noCSK_filt.RData")
#Make the data objects that will be needed
bg<-new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))
dataGMM<-new("GMMres", res=GMM, IDmap=GMM.ID, resFC=GMM.wFC)
#Choose the drug targets
targets.P<-list(cond1=c("MTOR_HUMAN"))
#Choose the drug treatments matching to the drug targets
#and match to what is present in the background network
data.P<-dataBycond(dataGMM, bg,scaled=TRUE,rowBycond=list(cond1=c("tor - ins")))
show(data.P)
#GMM results list, by condition: 603 peptides 
#Conditions: 1 
#Species perturbed: 
#  cond1.tor - ins : 34 
speciesP(data.P)
#$cond1
#$cond1$`tor - ins`
#[1] "PDCD4_HUMAN.S.76"   "4EBP2_HUMAN.Y.34"   "PDCD4_HUMAN.S.457"  "KIF4A_HUMAN.S.801" 
#[5] "NPM_HUMAN.S.70"     "KS6B1_HUMAN.S.447"  "TIF1B_HUMAN.S.594"  "TBCD4_HUMAN.T.752" 
#[9] "RPTOR_HUMAN.S.863"  "4EBP1_HUMAN.Y.34"   "PYR1_HUMAN.S.1859"  "TOIP1_HUMAN.S.156" 
#[13] "P53_HUMAN.S.314"    "AKT3_HUMAN.S.472"   "SRRM1_HUMAN.T.572"  "TIF1B_HUMAN.S.473" 
#[17] "NDRG1_HUMAN.S.332"  "TPIS_HUMAN.S.21"    "SYEP_HUMAN.S.882"   "HNRPK_HUMAN.S.379" 
#[21] "NDRG1_HUMAN.S.330"  "TBCD4_HUMAN.S.588"  "AKTS1_HUMAN.T.246"  "NPM_HUMAN.S.125"   
#[25] "SF01_HUMAN.S.80"    "BICD2_HUMAN.S.582"  "SCRIB_HUMAN.S.1486" "NU153_HUMAN.S.334" 
#[29] "TIF1B_HUMAN.S.50"   "TBCD4_HUMAN.S.591"  "NU107_HUMAN.S.11"   "TRIPC_HUMAN.S.991" 
#[33] "NOC2L_HUMAN.S.673"  "TIF1B_HUMAN.S.501" 
#Create the PKN list that will be used for optimisation
pknList<-buildNw(data.On=data.P, targets.On=targets.P, bg=bg,nK="no")
#[1] "Your data contains information about 34 sites, of which 15 are in your network"
#[1] "Your complete network contains 117 kinase/phosphatase substrate interactions (possibly non unique)"
#[1] "30 interactions to 15 integrator nodes were added"
show(pknList)
#Interactions: 78 
#Direct interactions with targets: 3 
#Species: 55 
#Sinks: 10 
#Integrators: 15 
#Intermediates: 30 
#To have a look at what your data contains: show(data.P)
#To have a look at the sites perturbed under each drug: speciesP(data.P)
#To have a look at your network: show(pknList)
#set the parameters for optimisation
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=5, nG=50)
save(file=paste("../cluster_scripts/torin/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=6, nG=50)
save(file=paste("../cluster_scripts/torin/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=7, nG=50)
save(file=paste("../cluster_scripts/torin/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=8, nG=50)
save(file=paste("../cluster_scripts/torin/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=3, nG=50)
save(file=paste("../cluster_scripts/torin/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))
optParam<-new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE,
              nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=4, nG=50)
save(file=paste("../cluster_scripts/torin/data4cluster_",resN(optParam),".RData", sep=""), 
     list=c("pknList","data.P", "targets.P", "optParam", "dataGMM"))

