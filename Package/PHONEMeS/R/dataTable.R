dataTable<-function(nodesP, drugs, dataGMM){
	#nodesP= character vector containing the identifiers of perturbed nodes
	#in format: UPID.res.site
	#drugs is a character vector that contains the drug treatments
	#that we want the data for (they must match to rownames of the GMM.res elements)
	#rsFC(dataGMM) is a list resulting from the GMM, with an element for each peptide,
	#and inside each element a data matrix with named rows corresponding to drugs
	#and named columns "clus", "status" and "FCvC" (at least) (containing respectively
	#the cluster from the GMM (P,C,I), 
	#the status (OK if p value of FC matches cluster assignment, 
		#FP if clus=P but p-val. not significant, 
		#FN if clus=C but p-val significant, I if clus=I) )
	#the FC vs control
	#
	#this checks that the format of nodesP is right, and excludes the entries 
	#where the format isn't right
	whichIDs<-grep("\\w+_\\w+.\\w{1}.\\d+",nodesP, perl=TRUE)
	if(!all(1:length(nodesP) %in% whichIDs)){
		cat("The following node identifiers are not in the right format (UPID.res.site) and will be excluded:")
		cat("\n")
		cat(nodesP[-whichIDs], sep=",")
	}
	nodesP<-nodesP[whichIDs]
	#this gets the data IDs for those sites
	nodesP.table<-cbind(as.character(IDmap(dataGMM)[which(IDmap(dataGMM)$S.cc %in% nodesP),"S.cc"]),as.character(IDmap(dataGMM)[which(IDmap(dataGMM)$S.cc %in% nodesP), "dataID"]))
	nodesP.table<-nodesP.table[!duplicated(nodesP.table),]
	#this prepares an empty data table
	nodesP.table<-cbind(nodesP.table, matrix(NA, nrow=dim(nodesP.table)[1], ncol=length(drugs)*3))
	colnames(nodesP.table)<-c("Name", "dataID", paste(rep(drugs, each=3), rep(c("clus", "status", "FC"), length(drugs)), sep="_"))
	#fill in the table
	for(i in 1:dim(nodesP.table)[1]){
		if(as.character(nodesP.table[i,2]) %in% names(resFC(dataGMM))){
			dataL<-resFC(dataGMM)[as.character(nodesP.table[i, "dataID"])]
			dataL<-dataL[!unlist(lapply(dataL,is.null))]
			colD<-c(which(colnames(dataL[[1]]) == "clus"), 
					which(colnames(dataL[[1]]) == "status"),
					ifelse("FCvC" %in% colnames(dataL[[1]]), which(colnames(dataL[[1]]) == "FCvC"), which(colnames(dataL[[1]]) == "FC")))
			nodesP.table[i,3:dim(nodesP.table)[2]]<-as.vector(t(dataL[[1]][drugs,colD]))
			}
	}
	notFound<-apply(nodesP.table, MARGIN=1, function(x){all(is.na(x[3:length(x)]))})
	nodesP.table<-nodesP.table[!notFound,]
	return(nodesP.table)
}
