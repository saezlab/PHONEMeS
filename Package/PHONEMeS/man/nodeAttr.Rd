\name{nodeAttr}
\alias{nodeAttr}
\title{
Writes node attributes for a network and data set
}
\description{
This function writes node attributes and data tables for a data set matched to a network
}
\usage{
nodeAttr(targets.On, nodesOF, drugsD, dataGMM, optParam, resFolder=getwd())
}
\arguments{
  \item{targets.On}{
a list containing IDs of K/P perturbed in each condition
}
  \item{drugsD}{
a character vector that contains the drug treatments that we want the data for (they must match to row names of the data frames in the \code{res(dataGMM)} element)
}
  \item{nodesOF}{
a list of nodes perturbed or not in a set of conditions, as produced by \code{nodesData}
}
  \item{optParam}{
an \code{OptParam} object
}
  \item{resFolder}{
a character vector, the path to the folder where results will be written
}
  \item{dataGMM}{
a \code{GMMres} object, with at least the elements \code{resFC} and \code{IDmap}. \code{resFC} is a list resulting from the GMM, with a named element for each peptide (where the names map to the "dataID" column in \code{IDmap(dataGMM)}, and inside each element a data matrix with named rows corresponding to drugs and named columns "clus", "status" and "FCvC", containing respectively the cluster from the GMM (P/C/I), the status (OK if p-value of fold change matches cluster assignment, FP if clus=P but p-val. not significant, FN if clus=C but p-val significant, I if clus=I), and the fold change vs control). \code{IDmap} is a data frame or matrix containing the map from site IDs to data IDs, with site IDs in a column named "S.cc" (mapping to nodesP), data ID in column named "dataID" (mapping to the names of \code{resFC}).
}
}
\value{
This function writes NA and DA files
}
\author{
C. Terfve
}
