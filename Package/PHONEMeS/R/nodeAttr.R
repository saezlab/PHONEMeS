nodeAttr<-function(targets.On, nodesOF, drugsD, dataGMM, optParam, resFolder=getwd()){
	#These write the nodes attributes file (1 general and 1 for each condition)
	write.table(rbind(c('Name', 'nodesP'),
	                  cbind(
              c(unlist(targets.On),unlist(nodesOF$On)), 
              c(rep("D",length(unlist(targets.On))), rep("P", length(unlist(nodesOF$On)))))),
            file=paste(resFolder,"AllNodes_nodesP_NA_p",resN(optParam),".txt", sep=""),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
	for(i in 1:length(targets.On)){
		write.table(cbind(c(targets.On[[i]],nodesOF$OnList[[i]], nodesOF$OffList[[i]], intersect(nodesOF$OnList[[i]], nodesOF$OffList[[i]])), 
			c(rep("D",length(targets.On[[i]])), rep("P", length(nodesOF$OnList[[i]])),rep("C", length(nodesOF$OffList[[i]])), 
			rep("B", length(intersect(nodesOF$OnList[[i]], nodesOF$OffList[[i]]))))),
            	file=paste(resFolder,"AllNodes_", i,"_NA_p",resN(optParam),".txt", sep=""),
            		sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    }
	##write the data table
	dT<-dataTable(nodesP=unlist(nodesOF$On), drugs=drugsD, dataGMM=dataGMM)
	write.table(dT, row.names=FALSE, col.names=TRUE, quote=FALSE, file=paste0(resFolder,"AllNodes_",resN(optParam),"_DA.txt"), sep="\t")
}