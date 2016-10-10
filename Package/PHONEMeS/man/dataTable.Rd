\name{dataTable}
\alias{dataTable}
\title{
Create a data table for sites in a network
}
\description{
This function creates a data table for a resulting network (i.e. the data that is used to annotate networks in terms of cluster (P/C), status (OK/FN/FP/I) and fold change sign, in Cytoscape)
}
\usage{
dataTable(nodesP, drugs, dataGMM)
}
\arguments{
  \item{nodesP}{
a character vector containing the substrate ID of data points as UPID.res.pos (this is to select sites that should be mapped)
}
  \item{drugs}{
a character vector that contains the drug treatments that we want the data for (they must match to row names of the data frames in the \code{res(dataGMM)} element)
}
  \item{dataGMM}{
a \code{GMMres} object, with at least the elements \code{resFC} and \code{IDmap}. \code{resFC} is a list resulting from the GMM, with a named element for each peptide (where the names map to the "dataID" column in \code{IDmap(dataGMM)}, and inside each element a data matrix with named rows corresponding to drugs and named columns "clus", "status" and "FCvC", containing respectively the cluster from the GMM (P/C/I), the status (OK if p-value of fold change matches cluster assignment, FP if clus=P but p-val. not significant, FN if clus=C but p-val significant, I if clus=I), and the fold change vs control). \code{IDmap} is a data frame or matrix containing the map from site IDs to data IDs, with site IDs in a column named "S.cc" (mapping to nodesP), data ID in column named "dataID" (mapping to the names of \code{resFC}).
}
}
\value{
This function returns a data matrix that contains, for each site and each drug used, the cluster, status and fold change information. This is meant to be used as a node attribute in Cytoscape, to visualise the results of an analysis.
}
\author{
C. Terfve
}
