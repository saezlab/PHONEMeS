\name{dataBycond}
\alias{dataBycond}
\title{
Function that re-shapes the data so that its entries match the targets list
}
\description{
Function that re-shapes the data so that its entries match the targets list. The targets list is a list with an element for each condition, where a condition is defined as a set of targets  that need to be on. This function creates a matching list that has the corresponding data for nodes that are reachable by the base network (before building of the custom network).
}
\usage{
dataBycond(dataGMM=dataGMM, bg=bg, scaled=TRUE, rowBycond=list(cond1=c("MTOR1 - Control2","MTOR2 - Control2")))
}
\arguments{
  \item{dataGMM}{
a \code{GMMres} object that contains at least the \code{res} and \code{IDmap}. \code{res} is the result of the Gaussian mixture modelling process, i.e. a list with an element for each peptide, each element containing a matrix with: the P/C cluster assignment ("clus"), log likelihood of P/C ("Indiv"), adjusted p value of FC drug vs ctrl ("FCvCaPval") and OK/FP/FN/I status ("status").
}
  \item{bg}{
a \code{KPSbg} object that contains all interactions that can be considered to make a background network, e.g. made from the \code{allD} data frame resulting from the database combination process. The original interactions data frame must have at least the following columns: the "K.ID" column (kinase Uniprot ID) and the "S.ID" column (substrate Uniprot ID), the "SID" column (unique character string identifiers for interactions), the "S.cc" column (substrateUPID.res.pos)
}
  \item{rowBycond}{
a list with an element for each condition (where a condition refers to a set of experimental treatments that have the same K/P targets in the networks), containing a vector of character strings that specify which rows of the \code{res(GMMres)} matrix should be used for this condition. 
}
  \item{scaled}{
logical: do we want the scores to be scaled between -1 and 1 by drug
}
}
\details{
This function reshapes the data in a \code{GMMres} object so as to have a list of lists where each element is a condition (set of K/P targets) and each sub element is a drug treatment acting on that set of targets.  If cond1 refers to a pair of MTOR inhibitors, then the targets list has an element cond1 that only contains \code{"MTOR_HUMAN"} and the data list created here has a slot called "cond1" that has 2 data frames, one for each of the 2 drugs that match this target.
}
\value{
This function returns a \code{GMMbyCond} object, which has 2 elements: a \code{dataBC} list with one element for each condition, and inside one element for each drug treatment that matches that condition; a \code{speciesP} element which is a list matching to the above, containing the nodes in the data (in the same format as S.cc, i.e. "UPID.res.pos") that are perturbed under the drug treatments used. The data in both of these elements only relates to those nodes that are reachable in the full background network.
}
\author{
C. Terfve
}
