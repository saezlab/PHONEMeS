data<-read.table("1199498s_TableS2.txt", sep="\t", header=TRUE)
dataIDs<-as.character(data$refseq_id)
data.sites<-as.character(data$abs_phosphosites)
data.IDmap<-read.csv("uniprot-yourlist%3AM201501042EL1V9A9QJ-2.csv", header=TRUE,
                     fill=TRUE)
data.IDtable<-cbind(
  as.character(data.IDmap$yourlist.M201501042EL1V9A9QJ[data.IDmap$Status == "reviewed"]), 
  as.character(data.IDmap$Entry.name[data.IDmap$Status == "reviewed"]))
#first I'll deconvolute run A and B, and the ctrl vs Ins vs Ins+Tor vs Ins+rapa
data<-cbind(data,paste(data$refseq_id, data$abs_phosphosites, sep="_"))
data.runA<-data[which(data$run == "A"),]
data.runB<-data[which(data$run == "B"),]
#A: 2475 rows
#B: 3531 rows
#now combine by peptide
data.2<-unique(as.character(data[,17]))
#3883
data.2<-data.frame(Scc=data.2, matrix(NA, nrow=length(data.2), ncol=8))
colnames(data.2)<-c("Scc", "refseq_id", "abs_phosphosites",
                    "ctrl_ins_log2.A","tor_ins_log2.A","rapa_ins_log2.A",
                    "ctrl_ins_log2.B","tor_ins_log2.B","rapa_ins_log2.B")
data.2$refseq_id<-as.character(data$refseq_id[match(data.2$Scc, data[,17])])
data.2$abs_phosphosites<-as.character(data$abs_phosphosites[match(data.2$Scc, data[,17])])
data.2$ctrl_ins_log2.A<-as.numeric(data.runA$ctrl_ins_log2[match(data.2$Scc,data.runA[,17])])
data.2$tor_ins_log2.A<-as.numeric(data.runA$tor_ins_log2[match(data.2$Scc,data.runA[,17])])
data.2$rapa_ins_log2.A<-as.numeric(data.runA$rapa_ins_log2[match(data.2$Scc,data.runA[,17])])
data.2$ctrl_ins_log2.B<-as.numeric(data.runB$ctrl_ins_log2[match(data.2$Scc,data.runB[,17])])
data.2$tor_ins_log2.B<-as.numeric(data.runB$tor_ins_log2[match(data.2$Scc,data.runB[,17])])
data.2$rapa_ins_log2.B<-as.numeric(data.runB$rapa_ins_log2[match(data.2$Scc,data.runB[,17])])
#
data.2.rep<-data.2[!is.na(data.2$rapa_ins_log2.A),]
data.2.rep<-data.2.rep[!is.na(data.2.rep$rapa_ins_log2.B),]
#excluding those peptides that do not have replicates I go from
#3883 peptides to 1609 peptides
plot(density(data.2.rep$tor_ins_log2.A, na.rm=TRUE), ylim=c(0,1.5), col="blue")
lines(density(data.2.rep$rapa_ins_log2.A, na.rm=TRUE), col="green")
lines(density(data.2.rep$ctrl_ins_log2.A, na.rm=TRUE))
#
plot(density(data.2.rep$tor_ins_log2.B, na.rm=TRUE), ylim=c(0,2), col="blue")
lines(density(data.2.rep$rapa_ins_log2.B, na.rm=TRUE), col="green")
lines(density(data.2.rep$ctrl_ins_log2.B, na.rm=TRUE))
#center/standardise
data.2.repN<-data.2.rep
data.2.repN$tor_ins_log2.A<-(data.2.repN$tor_ins_log2.A-mean(data.2.repN$tor_ins_log2.A))/sd(data.2.repN$tor_ins_log2.A)
data.2.repN$rapa_ins_log2.A<-(data.2.repN$rapa_ins_log2.A-mean(data.2.repN$rapa_ins_log2.A))/sd(data.2.repN$rapa_ins_log2.A)
data.2.repN$ctrl_ins_log2.A<-(data.2.repN$ctrl_ins_log2.A-mean(data.2.repN$ctrl_ins_log2.A))/sd(data.2.repN$ctrl_ins_log2.A)
data.2.repN$tor_ins_log2.B<-(data.2.repN$tor_ins_log2.B-mean(data.2.repN$tor_ins_log2.B))/sd(data.2.repN$tor_ins_log2.B)
data.2.repN$rapa_ins_log2.B<-(data.2.repN$rapa_ins_log2.B-mean(data.2.repN$rapa_ins_log2.B))/sd(data.2.repN$rapa_ins_log2.B)
data.2.repN$ctrl_ins_log2.B<-(data.2.repN$ctrl_ins_log2.B-mean(data.2.repN$ctrl_ins_log2.B))/sd(data.2.repN$ctrl_ins_log2.B)
#
plot(density(data.2.repN$tor_ins_log2.A, na.rm=TRUE), ylim=c(0,1.5), col="blue")
lines(density(data.2.repN$rapa_ins_log2.A, na.rm=TRUE), col="green")
lines(density(data.2.repN$ctrl_ins_log2.A, na.rm=TRUE))
plot(density(data.2.repN$tor_ins_log2.B, na.rm=TRUE), ylim=c(0,2), col="blue")
lines(density(data.2.repN$rapa_ins_log2.B, na.rm=TRUE), col="green")
lines(density(data.2.repN$ctrl_ins_log2.B, na.rm=TRUE))
#Average the replicates
tor_ins_log2<-data.2.repN[,c("Scc","refseq_id","abs_phosphosites")]
tor_ins_log2<-cbind(tor_ins_log2, rowMeans(data.2.repN[,c("tor_ins_log2.A","tor_ins_log2.B")]))
sum(tor_ins_log2[,4] > (median(tor_ins_log2[,4])+2.5*mad(tor_ins_log2[,4])))#31
sum(tor_ins_log2[,4] < (median(tor_ins_log2[,4])-2.5*mad(tor_ins_log2[,4])))#52
rapa_ins_log2<-data.2.repN[,c("Scc","refseq_id","abs_phosphosites")]
rapa_ins_log2<-cbind(rapa_ins_log2, rowMeans(data.2.repN[,c("rapa_ins_log2.A","rapa_ins_log2.B")]))
sum(rapa_ins_log2[,4] > (median(rapa_ins_log2[,4])+2.5*mad(rapa_ins_log2[,4])))#51
sum(rapa_ins_log2[,4] < (median(rapa_ins_log2[,4])-2.5*mad(rapa_ins_log2[,4])))#60
#plot the difference between the replicates divided by the mean, by peptides
plot(abs(data.2.repN[,"ctrl_ins_log2.A"]-data.2.repN[,"ctrl_ins_log2.B"])/(abs(rowMeans(data.2.repN[,c("ctrl_ins_log2.A", "ctrl_ins_log2.B")]))))
points(abs(data.2.repN[,"tor_ins_log2.A"]-data.2.repN[,"tor_ins_log2.B"])/(abs(rowMeans(data.2.repN[,c("tor_ins_log2.A", "tor_ins_log2.B")]))), col="blue")
points(abs(data.2.repN[,"rapa_ins_log2.A"]-data.2.repN[,"rapa_ins_log2.B"])/(abs(rowMeans(data.2.repN[,c("rapa_ins_log2.A", "rapa_ins_log2.B")]))), col="green")
dM<-cbind(abs(data.2.repN[,"ctrl_ins_log2.A"]-data.2.repN[,"ctrl_ins_log2.B"])/(abs(rowMeans(data.2.repN[,c("ctrl_ins_log2.A", "ctrl_ins_log2.B")]))),
          abs(data.2.repN[,"tor_ins_log2.A"]-data.2.repN[,"tor_ins_log2.B"])/(abs(rowMeans(data.2.repN[,c("tor_ins_log2.A", "tor_ins_log2.B")]))),
          abs(data.2.repN[,"rapa_ins_log2.A"]-data.2.repN[,"rapa_ins_log2.B"])/(abs(rowMeans(data.2.repN[,c("rapa_ins_log2.A", "rapa_ins_log2.B")]))))
#now plot the sum of these
plot(density(rowSums(dM)), xlim=c(0,100))
#I'll remove the above 10 
sum(rowSums(dM)>10)
#567
#I'll call P what's beyond 2.5 MAD
#I'll call OK what is either C or P AND with dM<10
sum(rowSums(dM)<10)#1042
#now I'll get a p-value for different from the mean
tor_ins_log2.p<-pnorm(tor_ins_log2[,4], mean=mean(tor_ins_log2[,4]), sd=sd(tor_ins_log2[,4]), lower.tail=(tor_ins_log2[,4] < mean(tor_ins_log2[,4])))
# if lower.tail was false then I need 1-p
tor_ins_log2.p[which(tor_ins_log2[,4] > mean(tor_ins_log2[,4]))]<-1-tor_ins_log2.p[which(tor_ins_log2[,4] > mean(tor_ins_log2[,4]))]
tor_ins_log2.ap<-p.adjust(tor_ins_log2.p)
sum(tor_ins_log2.p < 0.05)#108
sum(tor_ins_log2.ap < 0.05)#9
#get p/1-p  
tor_ins_log2.lo<-log(tor_ins_log2.p/(1-tor_ins_log2.p), base=10)
#still have the P cluster assignment as in the paper (2.5 MAD)
rapa_ins_log2.p<-pnorm(rapa_ins_log2[,4], mean=mean(rapa_ins_log2[,4]), sd=sd(rapa_ins_log2[,4]), lower.tail=(rapa_ins_log2[,4] < mean(rapa_ins_log2[,4])))
# if lower.tail was false then I need 1-p
rapa_ins_log2.p[which(rapa_ins_log2[,4] > mean(rapa_ins_log2[,4]))]<-1-rapa_ins_log2.p[which(rapa_ins_log2[,4] > mean(rapa_ins_log2[,4]))]
rapa_ins_log2.ap<-p.adjust(rapa_ins_log2.p)
sum(rapa_ins_log2.p < 0.05)#115
sum(rapa_ins_log2.ap < 0.05)#8
#get p/1-p  
rapa_ins_log2.lo<-log(rapa_ins_log2.p/(1-rapa_ins_log2.p), base=10)
hist(rapa_ins_log2.p)
hist(rapa_ins_log2.lo)
rapa_ins_log2<-cbind(rapa_ins_log2, rapa_ins_log2.p, rapa_ins_log2.ap, rapa_ins_log2.lo)
tor_ins_log2<-cbind(tor_ins_log2, tor_ins_log2.p, tor_ins_log2.ap, tor_ins_log2.lo)
rownames(rapa_ins_log2)<-rapa_ins_log2[,"Scc"]
rownames(tor_ins_log2)<-tor_ins_log2[,"Scc"]
tor_ins_log2.c<-rep("C", dim(tor_ins_log2)[1])
tor_ins_log2.c[which(tor_ins_log2[,4] > (median(tor_ins_log2[,4])+2.5*mad(tor_ins_log2[,4])))]<-"P"
tor_ins_log2.c[which(tor_ins_log2[,4] < (median(tor_ins_log2[,4])-2.5*mad(tor_ins_log2[,4])))]<-"P"
table(tor_ins_log2.c)
#C    P 
#1526   83 
rapa_ins_log2.c<-rep("C", dim(rapa_ins_log2)[1])
rapa_ins_log2.c[which(rapa_ins_log2[,4] > (median(rapa_ins_log2[,4])+2.5*mad(rapa_ins_log2[,4])))]<-"P"
rapa_ins_log2.c[which(rapa_ins_log2[,4] < (median(rapa_ins_log2[,4])-2.5*mad(rapa_ins_log2[,4])))]<-"P"
table(rapa_ins_log2.c)
#C    P 
#1498  111
#status
tor_ins_log2.s<-rep("OK", dim(tor_ins_log2)[1])
rapa_ins_log2.s<-rep("OK", dim(rapa_ins_log2)[1])
names(tor_ins_log2.s)<-tor_ins_log2[,"Scc"]
names(rapa_ins_log2.s)<-rapa_ins_log2[,"Scc"]
tor_ins_log2.s[as.character(data.2.repN[which(rowSums(dM)>10), "Scc"])]<-"FP"
rapa_ins_log2.s[as.character(data.2.repN[which(rowSums(dM)>10), "Scc"])]<-"FP"
table(tor_ins_log2.s)
#FP   OK 
#567 1042 
#
rapa_ins_log2<-cbind(rapa_ins_log2, rapa_ins_log2.c, rapa_ins_log2.s)
tor_ins_log2<-cbind(tor_ins_log2, tor_ins_log2.c, tor_ins_log2.s)
for(i in 1:dim(rapa_ins_log2)[1]){
  if(rapa_ins_log2[i,"rapa_ins_log2.c"] == "C") rapa_ins_log2[i,"rapa_ins_log2.s"]<-"OK"
  if(tor_ins_log2[i,"tor_ins_log2.c"] == "C") tor_ins_log2[i,"tor_ins_log2.s"]<-"OK"
}
#
GMM<-vector("list", dim(tor_ins_log2)[1])
names(GMM)<-as.character(tor_ins_log2[,"Scc"])
for(i in 1:length(GMM)){
  GMM[[i]]<-rbind(
    c(as.character(rapa_ins_log2[names(GMM)[i],"rapa_ins_log2.lo"]), 
      as.character(rapa_ins_log2[names(GMM)[i],"rapa_ins_log2.c"]),
      as.character(rapa_ins_log2[names(GMM)[i],"rapa_ins_log2.ap"]),
      as.character(rapa_ins_log2[names(GMM)[i],"rapa_ins_log2.s"])),
    c(as.character(tor_ins_log2[names(GMM)[i],"tor_ins_log2.lo"]), 
      as.character(tor_ins_log2[names(GMM)[i],"tor_ins_log2.c"]),
      as.character(tor_ins_log2[names(GMM)[i],"tor_ins_log2.ap"]),
      as.character(tor_ins_log2[names(GMM)[i],"tor_ins_log2.s"]))
  )
  colnames(GMM[[i]])<-c("Indiv", "clus","FCvCaPval","status")
  rownames(GMM[[i]])<-c("rapa - ins", "tor - ins")
}
#Remove unmapped
GMM.ID.2<-cbind(as.character(tor_ins_log2[,"Scc"]),
              as.character(data.IDtable[match(tor_ins_log2[,"refseq_id"], as.character(data.IDtable[,1])),2]),
              as.character(tor_ins_log2[,"abs_phosphosites"])
)
GMM<-GMM[GMM.ID.2[!is.na(GMM.ID.2[,2]),1]]
#
GMM.wFC<-vector("list", dim(tor_ins_log2)[1])
names(GMM.wFC)<-as.character(tor_ins_log2[,"Scc"])
for(i in 1:length(GMM.wFC)){
  GMM.wFC[[i]]<-rbind(
    c(as.character(rapa_ins_log2[names(GMM.wFC)[i],"rapa_ins_log2.lo"]), 
      as.character(rapa_ins_log2[names(GMM.wFC)[i],"rapa_ins_log2.c"]),
      as.character(rapa_ins_log2[names(GMM.wFC)[i],"rapa_ins_log2.ap"]),
      as.character(rapa_ins_log2[names(GMM.wFC)[i],"rapa_ins_log2.s"]),
      as.character(rapa_ins_log2[names(GMM.wFC)[i],4])),
    c(as.character(tor_ins_log2[names(GMM.wFC)[i],"tor_ins_log2.lo"]), 
      as.character(tor_ins_log2[names(GMM.wFC)[i],"tor_ins_log2.c"]),
      as.character(tor_ins_log2[names(GMM.wFC)[i],"tor_ins_log2.ap"]),
      as.character(tor_ins_log2[names(GMM.wFC)[i],"tor_ins_log2.s"]),
      as.character(tor_ins_log2[names(GMM.wFC)[i],4]))
  )
  colnames(GMM.wFC[[i]])<-c("Indiv", "clus","FCvCaPval","status", "FCvC")
  rownames(GMM.wFC[[i]])<-c("rapa - ins", "tor - ins")
}
GMM.wFC<-GMM.wFC[GMM.ID.2[!is.na(GMM.ID.2[,2]),1]]
#
rm(GMM.ID.2)
GMM.ID<-data.frame(dataID=as.character(GMM.ID[,1]), 
                   UPID=as.character(GMM.ID[,2]),
                   site=as.character(GMM.ID[,3]),
                   pos=as.numeric(GMM.ID[,4]),
                   res=as.character(GMM.ID[,5]),
                   S.cc=as.character(GMM.ID[,6]))
#
save(list=c("GMM.ID", "GMM","GMM.wFC", "rapa_ins_log2", "tor_ins_log2", "data.2.repN"), file="data_Hsu_PHONEMeS_Feb13.RData")
