
# Set working directory to directory of this script (in RStudio)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# If you are using this file as Source, use:
# setwd(getSrcDirectory()[1])



load("../data/MSdata.RData")
library(limma)
###########Normalise and apply the linear model
d12.log<-log(d12, base=10)
d12.log.dN<-normalizeBetweenArrays(as.matrix(d12.log), method="quantile")
rownames(d12.log.dN)<-rownames(d12)
colnames(d12.log.dN)<-colnames(d12)
conditions1<-c("AKT1","AKT2","CAMK1","CAMK2","EGFR1" ,"EGFR2",
               "ERK1","ERK2" ,"MEK1","MEK2","MTOR1","MTOR2","Control1")
conditions2<-c("P70S6K1","P70S6K2","PI3K1" ,"PI3K2","PKC1","PKC2",
               "PP2A1","PP2A2","ROCK1","ROCK2","Control2")
design6<-matrix(0,ncol=27, nrow=144)
colnames(design6)<-c("Control",conditions1[1:12], conditions2[1:10], "E1", "E1a", "E2", "E2a")
rownames(design6)<-colnames(d12.log.dN)
design6[c("Control1","Control1.1","Control1.2",
          "Control1a","Control1a.1","Control1a.2",
          "Control2","Control2.1","Control2.2",
          "Control2a","Control2a.1","Control2a.2"),"Control"]<-1
design6[c("AKT1","AKT1.1","AKT1.2","AKT1a","AKT1a.1","AKT1a.2"),"AKT1"]<-1
design6[c("AKT2","AKT2.1","AKT2.2",
          "AKT2a","AKT2a.1","AKT2a.2"),"AKT2"]<-1
design6[c("CAMK1","CAMK1.1","CAMK1.2",
          "CAMK1a","CAMK1a.1","CAMK1a.2"),"CAMK1"]<-1
design6[c( "CAMK2","CAMK2.1","CAMK2.2",
           "CAMK2a","CAMK2a.1","CAMK2a.2"),"CAMK2"]<-1
design6[c("EGFR1","EGFR1.1","EGFR1.2",
          "EGFR1a","EGFR1a.1","EGFR1a.2"),"EGFR1"]<-1
design6[c("EGFR2","EGFR2.1","EGFR2.2",
          "EGFR2a","EGFR2a.1","EGFR2a.2"),"EGFR2"]<-1
design6[c("ERK1","ERK1.1","ERK1.2",
          "ERK1a","ERK1a.1","ERK1a.2"),"ERK1"]<-1 
design6[c("ERK2","ERK2.1","ERK2.2" ,
          "ERK2a","ERK2a.1","ERK2a.2"),"ERK2"]<-1
design6[c("MEK1","MEK1.1","MEK1.2",
          "MEK1a","MEK1a.1","MEK1a.2" ),"MEK1"]<-1
design6[c("MEK2","MEK2.1","MEK2.2",
          "MEK2a","MEK2a.1","MEK2a.2"),"MEK2"]<-1
design6[c("MTOR1","MTOR1.1","MTOR1.2",
          "MTOR1a","MTOR1a.1","MTOR1a.2"),"MTOR1"]<-1
design6[c("MTOR2","MTOR2.1","MTOR2.2",
          "MTOR2a","MTOR2a.1","MTOR2a.2"),"MTOR2"]<-1
design6[c("P70S6K1","P70S6K1.1","P70S6K1.2",
          "P70S6K1a","P70S6K1a.1","P70S6K1a.2"),"P70S6K1"]<-1
design6[c("P70S6K2","P70S6K2.1","P70S6K2.2" ,
          "P70S6K2a","P70S6K2a.1","P70S6K2a.2"),"P70S6K2"]<-1
design6[c("PI3K1","PI3K1.1","PI3K1.2",
          "PI3K1a","PI3K1a.1","PI3K1a.2"),"PI3K1"]<-1
design6[c("PI3K2","PI3K2.1","PI3K2.2" ,
          "PI3K2a","PI3K2a.1","PI3K2a.2"),"PI3K2"]<-1
design6[c("PKC1","PKC1.1","PKC1.2",
          "PKC1a","PKC1a.1","PKC1a.2"),"PKC1"]<-1
design6[c("PKC2","PKC2.1","PKC2.2",
          "PKC2a","PKC2a.1","PKC2a.2"),"PKC2"]<-1
design6[c("PP2A1","PP2A1.1","PP2A1.2",
          "PP2A1a","PP2A1a.1","PP2A1a.2"),"PP2A1"]<-1
design6[c("PP2A2","PP2A2.1","PP2A2.2",
          "PP2A2a","PP2A2a.1","PP2A2a.2"),"PP2A2"]<-1
design6[c("ROCK1","ROCK1.1","ROCK1.2",
          "ROCK1a","ROCK1a.1","ROCK1a.2"),"ROCK1"]<-1
design6[c("ROCK2","ROCK2.1","ROCK2.2",
          "ROCK2a","ROCK2a.1","ROCK2a.2"),"ROCK2"]<-1
design6[c(1:39),"E1"]<-1
design6[c(40:78), "E1a"]<-1
design6[c(79:111),"E2"]<-1
design6[c(112:144),"E2a"]<-1
#
dC<-duplicateCorrelation(d12.log.dN, design=design6, block=rep(1:48, each=3))
fit6.b<-lmFit(d12.log.dN, design6, block=rep(1:48, each=3), correlation=dC$consensus.correlation)
#Coefficients not estimable: E2a 
#Partial NA coefficients for 12240 probe(s) 
dim(d12.log.dN)#[1] 12266   144
dim(fit6.b$coefficients)#[1] 12266    27
cm6<-makeContrasts(AKT1-Control,
                   AKT2-Control,
                   CAMK1-Control,
                   CAMK2-Control,
                   EGFR1-Control,
                   EGFR2-Control,
                   ERK1-Control,
                   ERK2-Control,
                   MEK1-Control,
                   MEK2-Control,
                   MTOR1-Control,
                   MTOR2-Control,
                   P70S6K1-Control,
                   P70S6K2-Control,
                   PI3K1-Control,
                   PI3K2-Control,
                   PKC1-Control,
                   PKC2-Control,
                   PP2A1-Control,
                   PP2A2-Control,
                   ROCK1-Control,
                   ROCK2-Control
                   , levels=design6)
fit6.b.2<-contrasts.fit(fit6.b, cm6)
fit6.b.2<-eBayes(fit6.b.2)
results6.b<-decideTests(fit6.b.2)
summary(results6.b)
#   AKT1 - Control AKT2 - Control CAMK1 - Control CAMK2 - Control
#-1           1646            628             321             203
#0            7783           8817            9246            9317
#1             232            216              94             141
#   EGFR1 - Control EGFR2 - Control ERK1 - Control ERK2 - Control
#-1             338             135            309            241
#0             9288            9501           9279           9384
#1               35              25             73             36
#   MEK1 - Control MEK2 - Control MTOR1 - Control MTOR2 - Control
#-1             67             33             281             319
#0            9582           9614            9321            9235
#1              12             14              59             107
#   P70S6K1 - Control P70S6K2 - Control PI3K1 - Control PI3K2 - Control
#-1               342               311             289             371
#0               9221              9210            9239            9126
#1                 98               140             133             164
#   PKC1 - Control PKC2 - Control PP2A1 - Control PP2A2 - Control
#-1            243            204            2278             435
#0            9357           9412            5997            8995
#1              61             45            1386             231
#   ROCK1 - Control ROCK2 - Control
#-1             238             103
#0             9399            9533
#1               24              25
cm6.dvd<-makeContrasts(AKT1-AKT2,
                       CAMK1-CAMK2,
                       EGFR1-EGFR2,
                       ERK1-ERK2,
                       MEK1-MEK2,
                       MTOR1-MTOR2,
                       P70S6K1-P70S6K2,
                       PI3K1-PI3K2,
                       PKC1-PKC2,
                       PP2A1-PP2A2,
                       ROCK1-ROCK2
                       , levels=design6)
fit6.b.2.dvd<-contrasts.fit(fit6.b, cm6.dvd)
fit6.b.2.dvd<-eBayes(fit6.b.2.dvd)
results6.b.dvd<-decideTests(fit6.b.2.dvd)
summary(results6.b.dvd)
#   AKT1 - AKT2 CAMK1 - CAMK2 EGFR1 - EGFR2 ERK1 - ERK2 MEK1 - MEK2
#-1         159            39            33         138           4
#0         9482          9614          9600        9246        9655
#1           20             8            28         277           2
#   MTOR1 - MTOR2 P70S6K1 - P70S6K2 PI3K1 - PI3K2 PKC1 - PKC2 PP2A1 - PP2A2
#-1            20               112            16         114          1900
#0           9615              9435          9614        9498          6466
#1             26               114            31          49          1295
#   ROCK1 - ROCK2
#-1           105
#0           9514
#1             42
limma.FCvCtrl.d6<-fit6.b.2$coefficients
limma.FCvCtrl.d6.pval<-fit6.b.2$p.value
limma.FCvCtrl.d6.Adjpval<-apply(limma.FCvCtrl.d6.pval, MARGIN=2, function(x){p.adjust(x,method="BH")})
NAs<-apply(limma.FCvCtrl.d6, MARGIN=1, function(x){all(is.na(x))})
limma.FCvCtrl.d6.NAr<-limma.FCvCtrl.d6[!NAs,]
limma.FCvCtrl.d6.pval.NAr<-limma.FCvCtrl.d6.pval[!NAs,]
limma.FCvCtrl.d6.Adjpval.NAr<-limma.FCvCtrl.d6.Adjpval[!NAs,]
#
limma.alone.d6<-fit6.b$coefficients
limma.alone.d6<-limma.alone.d6[,c(1:23)]
NAs<-apply(limma.alone.d6, MARGIN=1, function(x){all(is.na(x))})
limma.alone.d6.NAr<-limma.alone.d6[!NAs,]
#remove PP2A1 and 2 (these 2 drugs weren't used in the GMM and analysis
limma.alone.d6.NAr<-limma.alone.d6.NAr[,-match(c("PP2A1","PP2A2"), colnames(limma.alone.d6.NAr))]
limma.FCvCtrl.d6.NAr<-limma.FCvCtrl.d6.NAr[,-match(c("PP2A1 - Control","PP2A2 - Control"),colnames(limma.FCvCtrl.d6.NAr))]
limma.FCvCtrl.d6.pval.NAr<-limma.FCvCtrl.d6.pval.NAr[,-match(c("PP2A1 - Control","PP2A2 - Control"),colnames(limma.FCvCtrl.d6.pval.NAr))]
limma.FCvCtrl.d6.Adjpval.NAr<-limma.FCvCtrl.d6.Adjpval.NAr[,-match(c("PP2A1 - Control","PP2A2 - Control"),colnames(limma.FCvCtrl.d6.Adjpval.NAr))]
###########Gaussian Mixture Modelling
#Note that each run of the GMM will produce slightly different results because 
#of the stochastic optimisation in the GMM fitting
library(mclust)
#
gm.L.6<-rep(NA, dim(limma.alone.d6.NAr)[1])
gm.nobs<-rep(NA, dim(limma.alone.d6.NAr)[1])
gm.logL.L.6<-rep(NA, dim(limma.alone.d6.NAr)[1])
gm.bic.L.6<-rep(NA, dim(limma.alone.d6.NAr)[1])
gm.mu.L.6<-list(NA)
gm.var.L.6<-list(NA)
gm.proba.L.6<-list(NA)
gm.class.L.6<-matrix(NA, ncol=21, nrow=dim(limma.alone.d6.NAr)[1])
colnames(gm.class.L.6)<-colnames(limma.alone.d6.NAr)
rownames(gm.class.L.6)<-rownames(limma.alone.d6.NAr)
for(i in 1:dim(limma.alone.d6.NAr)[1]){
  gm.nobs[i]<-sum(is.na(limma.alone.d6.NAr[i,]))
  if(gm.nobs[i] < 11){
    gm<-Mclust(limma.alone.d6.NAr[i,!is.na(limma.alone.d6.NAr[i,])])
    gm.L.6[i]<-gm$G
    gm.logL.L.6[i]<-gm$loglik
    gm.bic.L.6[i]<-gm$bic
    gm.mu.L.6[[i]]<-gm$parameters$mean
    gm.var.L.6[[i]]<-gm$parameters$variance
    gm.proba.L.6[[i]]<-gm$z
    gm.class.L.6[i,names(gm$classification)]<-gm$classification
  }  
}
#these are the number of peptides that are best fitted with 1 to 9, components
table(gm.L.6)
#gm.L.6
#   1    2    3    4    5    6    7    8    9 
#7263 3183  670  230  132   67   50   34   25
k2.6<-which(gm.L.6 == 2)
#remove those peptides where the control is NA (because I can't assign a control cluster)
ctrl.k2.6<-limma.alone.d6.NAr[k2.6,1]
length(k2.6)#3183
length(which(ctrl.k2.6 != 0))#2784
k2.6<-k2.6[which(ctrl.k2.6 != 0)]
k2.6<-k2.6[which(rownames(limma.alone.d6.NAr)[k2.6] %in% rownames(limma.FCvCtrl.d6.Adjpval.NAr))]
length(k2.6)#2627
for(i in 1:length(k2.6)){
  if(length(gm.var.L.6[[k2.6[i]]]$sigmasq) == 1){
    gm.var.L.6[[k2.6[i]]]$sigmasq<-rep(gm.var.L.6[[k2.6[i]]]$sigmasq, 2)
  }
}
#This function checks for the overlap of distributions (to avoid largely overlapping C/P clusters)
min(limma.alone.d6.NAr[k2.6,], na.rm=TRUE)#[1] -8.125338
max(limma.alone.d6.NAr[k2.6,], na.rm=TRUE)#[1] -1.628209
OC<-rep(NA, length(k2.6))
xs <- seq(-10,+3, .01)
min.f1f2 <- function(x, mu1, mu2, sd1, sd2) {
  f1 <- dnorm(x, mean=mu1, sd=sd1)
  f2 <- dnorm(x, mean=mu2, sd=sd2)
  pmin(f1, f2)
}
for(i in 1:length(k2.6)){
  OC[i] <- sum(min.f1f2(xs, mu1=gm.mu.L.6[[k2.6[i]]][1], mu2=gm.mu.L.6[[k2.6[i]]][2], sd1=sqrt(gm.var.L.6[[k2.6[i]]]$sigmasq[1]), sd2=sqrt(gm.var.L.6[[k2.6[i]]]$sigmasq[2])))
}
length(OC)#2627
length(which(OC < 10))#2376
#compute the p value, for each obs, of belonging to each distrib.
pCalc.list<-list(NA)
for(i in 1:length(k2.6)){
  dataVec<-limma.alone.d6.NAr[rownames(gm.class.L.6)[k2.6[i]],]
  p1<-rep(NA, length(dataVec))
  p2<-rep(NA, length(dataVec))
  for(o in 1:length(dataVec)){
    if(!is.na(dataVec[o])){
      if(dataVec[o] < gm.mu.L.6[[k2.6[i]]][1]){
        p1[o]<-pnorm(dataVec[o], mean=gm.mu.L.6[[k2.6[i]]][1], sd=sqrt(gm.var.L.6[[k2.6[i]]]$sigmasq[1]), lower.tail=TRUE)
      }else{
        p1[o]<-pnorm(dataVec[o], mean=gm.mu.L.6[[k2.6[i]]][1], sd=sqrt(gm.var.L.6[[k2.6[i]]]$sigmasq[1]), lower.tail=FALSE)    
      }
      if(dataVec[o] < gm.mu.L.6[[k2.6[i]]][2]){
        p2[o]<-pnorm(dataVec[o], mean=gm.mu.L.6[[k2.6[i]]][2], sd=sqrt(gm.var.L.6[[k2.6[i]]]$sigmasq[2]), lower.tail=TRUE)  
      }else{
        p2[o]<-pnorm(dataVec[o], mean=gm.mu.L.6[[k2.6[i]]][2], sd=sqrt(gm.var.L.6[[k2.6[i]]]$sigmasq[2]), lower.tail=FALSE)      
      }
    }
    pCalc.list[[i]]<-cbind(p1,p2) 
  }
}
names(pCalc.list)<-rownames(gm.class.L.6)[k2.6]
#compute the cluster, simply based on which p value is the smaller
pCalc.list.clus<-list(NA)
for(i in 1:length(pCalc.list)){
  mat<-apply(pCalc.list[[i]],MARGIN=1, function(x){order(x)[1]})
  pCalc.list.clus[[i]]<-mat
}
names(pCalc.list.clus)<-rownames(gm.class.L.6)[k2.6]
for(i in 1:length(pCalc.list.clus)){
  part1<-pCalc.list.clus[[i]]
  part1[which(part1 == part1[1])]<-"C"
  part1[which(part1 != part1[1])]<-"P"
  pCalc.list.clus[[i]]<-part1
}
length(pCalc.list.clus)#[1] 2254
pCalc.list.clus.p<-pCalc.list.clus
for(i in 1:length(pCalc.list.clus)){
  p<-limma.FCvCtrl.d6.Adjpval.NAr[names(pCalc.list.clus)[i],]
  pCalc.list.clus.p[[i]]<-cbind(pCalc.list.clus.p[[i]], c(NA, p))
}
#This function checks for the classification accuracy of the GMM, based on agreement
#with the adjusted p-values from limma, for different values of the boundaries for 
#the log ratios (boundary) and different amounts of maximum overlap of the density curves
#of the 2 clusters (whichC)
classAcc5<-function(boundary=0.5, lLmax=+Inf, whichC=NA){
  pCalc.list.temp<-pCalc.list
  pCalc.list.clus.p.temp<-pCalc.list.clus.p
  logL<-unlist(gm.logL.L.6[k2.6])
  if(!is.na(whichC[1])){
    pCalc.list.temp<-pCalc.list.temp[whichC]
    pCalc.list.clus.p.temp<-pCalc.list.clus.p.temp[whichC]
    logL<-unlist(gm.logL.L.6[k2.6[whichC]])
  }
  pCalc.list.clusRatio.2<-pCalc.list.temp
  for(i in 1:length(pCalc.list.temp)){
    Corder2<-order(pCalc.list.temp[[i]][1,], decreasing=TRUE)
    pCalc.list.clusRatio.2[[i]][2:21]<-log(pCalc.list.temp[[i]][2:21,Corder2[1]]/pCalc.list.temp[[i]][2:21,Corder2[2]], base=10)
    clus<-rep("I", 21)
    clus[which(pCalc.list.clusRatio.2[[i]] > boundary)]<-"C"
    clus[which(pCalc.list.clusRatio.2[[i]] < -boundary)]<-"P"
    pCalc.list.clusRatio.2[[i]]<-cbind(pCalc.list.clusRatio.2[[i]][,1], clus)
  }
  table(unlist(lapply(pCalc.list.clusRatio.2, function(x){return(x[,2])})))
  for(i in 1:length(pCalc.list.clusRatio.2)){
    status<-rep("OK", 21)
    pCalc.list.clusRatio.2[[i]]<-cbind(pCalc.list.clusRatio.2[[i]], pCalc.list.clus.p.temp[[i]][,2])
    for(j in 1:dim(pCalc.list.clusRatio.2[[i]])[1]){
      if(any(is.na(pCalc.list.clusRatio.2[[i]][j,]))){
        status[j]<-NA
      }else{
        if(pCalc.list.clusRatio.2[[i]][j,2] == "P" && as.numeric(pCalc.list.clusRatio.2[[i]][j,3]) > 0.05){
          status[j]<-"FP"
        }
        if(pCalc.list.clusRatio.2[[i]][j,2] == "C" && as.numeric(pCalc.list.clusRatio.2[[i]][j,3]) < 0.05){
          status[j]<-"FN"
        }
        if(pCalc.list.clusRatio.2[[i]][j,2] == "I"){
          status[j]<-"I"
        }
      }
    }
    pCalc.list.clusRatio.2[[i]]<-cbind(pCalc.list.clusRatio.2[[i]], status) 
  }
  classifAcc.calc.I.2<-lapply(pCalc.list.clusRatio.2, function(x){return(x[,4])})
  classifAcc.calc.I.2<-classifAcc.calc.I.2[which(logL < lLmax)]
  print(length(classifAcc.calc.I.2))
  print(table(unlist(classifAcc.calc.I.2)))
  return(pCalc.list.clusRatio.2)
}
test2<-classAcc5(boundary=0.5, whichC=which(OC < 10))
#[1] [1] 2376
#   FN    FP     I    OK 
#  338  8725  1102 37355 
GMM.res<-test2
table(unlist(lapply(GMM.res, function(x){length(intersect(which(x[,"clus"] == "P"), which(x[,"status"] == "OK")))})))
#   0    1    2    3    4    5    6    7    8    9   10   11   12 
#1269  497  250  108   68   47   43   36   16   27    8    1    6  
#this might not be too bad
#497 sites are perturbed under 1 condition
##########IDs

########
#Fixing for the missing data

GMM.ID <- as.data.frame(matrix(, nrow = length(GMM.res), ncol = 6))
colnames(GMM.ID) <- c("dataID", "UPID", "site", "pos", "res", "S.cc")
GMM.ID$dataID <- names(GMM.res)
for(i in 1:nrow(GMM.ID)){
  geneID <- strsplit(GMM.ID$dataID[i], split = " ")[[1]][1]
  siteID <- strsplit(strsplit(GMM.ID$dataID[i], split = " ")[[1]][2], split = "-")[[1]][2]
  resID <- strsplit(siteID, split = "")[[1]][1]
  posID <- ""
  for(j in 2:length(strsplit(siteID, split = "")[[1]])){
    posID <- paste(posID, strsplit(siteID, split = "")[[1]][j], sep = "")
  }
  
  GMM.ID$UPID[i] <- paste(geneID, "_HUMAN", sep = "")
  GMM.ID$site[i] <- siteID
  GMM.ID$pos[i] <- posID
  GMM.ID$res[i] <- resID
  GMM.ID$S.cc[i] <- paste(GMM.ID$UPID[i], ".", resID, ".", posID, sep = "")
}

GMM.res.no.FC <- GMM.res
GMM.res <- list()
for(i in 1:length(GMM.res.no.FC)){
  tempGMM <- matrix(, nrow = dim(GMM.res.no.FC[[i]])[1], ncol = dim(GMM.res.no.FC[[i]])[2]+1)
  rownames(tempGMM) <- rownames(GMM.res.no.FC[[i]])
  tempGMM[, 1:dim(tempGMM)[2]-1] <- GMM.res.no.FC[[i]]
  rows <- rownames(GMM.res.no.FC[[i]])
  for(j in 2:nrow(tempGMM)){
    colID <- which(colnames(limma.FCvCtrl.d6.NAr)==rows[j])
    if(!is.null(colID)){
      tempGMM[j, dim(tempGMM)[2]] <- limma.FCvCtrl.d6.NAr[i, colID]
    }
    else{
      tempGMM[j, dim(tempGMM)[2]] <- NA
    }
  }
  colnames(tempGMM) <- c(colnames(GMM.res.no.FC[[i]]), "FCvC")
  GMM.res[[length(GMM.res)+1]] <- tempGMM
}

names(GMM.res) <- names(GMM.res.no.FC)

#saving the main results
save(list=c("GMM.res", "GMM.res.noFC", "GMM.res.ID", "d12", "d12.log.dN", "limma.FCvCtrl.d6.NAr", "limma.FCvCtrl.d6.pval.NAr", "limma.alone.d6.NAr"),file="../data/dataObjects_2_2.RData")


# FROM HERE ON NO NEED TO CONTINUE EXECUTION
# FROM HERE ON, THE CODE DOES NOT EXECUTE WITHOUT ERROR:
# object 'data2' not found

GMM.res.ID.2<-data.frame(dataID=as.character(data2$psite), UPID=as.character(data2$acc_no))
GMM.res.ID.2$site<-rep(NA, length(data2$psite))
DataSite<-as.character(data2$psite)
DataSite<-sub(DataSite, pattern="\\w+\\W", replacement="",perl=TRUE)
DataSite<-sub(DataSite,pattern="\\(z= \\d+\\)", replacement="", perl=TRUE)
DataSite<-gsub(DataSite,pattern="\\W", replacement="", perl=TRUE)
DataSite<-gsub(DataSite,pattern="Oxi", replacement="", perl=TRUE)
DataSite<-strsplit(DataSite, split="p")
DataSite<-lapply(DataSite, function(x){return(x[2:length(x)])})
names(DataSite)<-as.character(data2$psite)
for(i in 1:length(DataSite)){
  GMM.res.ID.2$site[i]<-DataSite[[i]][1]
  if(length(DataSite[[i]]) != 1){
    for(n in 2:length(DataSite[[i]])){
      rows2paste<-GMM.res.ID.2[i,]
      rows2paste$site<-DataSite[[i]][n]
      GMM.res.ID.2<-rbind(GMM.res.ID.2, rows2paste)
    }
  }  
}
GMM.res.ID.2$res<-gsub(GMM.res.ID.2$site, pattern="\\d+", replacement="", perl=TRUE)
GMM.res.ID.2$pos<-gsub(GMM.res.ID.2$site, pattern="\\D", replacement="", perl=TRUE)
GMM.res.ID.2$cc<-paste(GMM.res.ID.2$UPID, GMM.res.ID.2$res, GMM.res.ID.2$pos, sep=".")
length(unique(GMM.res.ID.2$cc))#[1] 11629
###
#this is all that's needed for the PHONEMeS part
save(list=c("GMM.res", "GMM.res.ID.2"),file="../data/GMM_noFC.RData") # This would explain the missing .RData file in the other script
#this more complete set of objects also has the fold changes
for(i in 1:length(GMM.res)){
  GMM.res[[i]]<-cbind(GMM.res[[i]], 
                      c(limma.FCvCtrl.d6.NAr[names(GMM.res)[i], match(rownames(GMM.res[[i]]), colnames(limma.FCvCtrl.d6.NAr))]))
  colnames(GMM.res[[i]])[c(1,3,5)]<-c("Indiv", "FCvCaPval", "FCvC")
}
save(list=c("GMM.res", "GMM.res.ID.2", "d12", "d12.log.dN", "limma.FCvCtrl.d6.NAr", "limma.FCvCtrl.d6.pval.NAr", "limma.alone.d6.NAr"),file="../data/dataObjects_2_2.RData")
####
