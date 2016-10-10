\name{buildNw}
\alias{buildNw}
\title{
Make a starting network object
}
\description{
This function makes a starting network object, it just calls the \code{buildPKN} and \code{PKNlist} functions
}
\usage{
buildNw(data.On,targets.On, bg,nK=c("all","no", "drugs2data", "data"))
}
\arguments{
  \item{data.On}{
a \code{GMMbyCond} object as created by \code{dataBycond}
}
  \item{targets.On}{
a list of drug targets (Uniprot identifiers matching to the "K.ID" column in the background network data frame)	
}
  \item{bg}{
a \code{KPSbg} object that contains all interactions that can be considered to make a background network, e.g. made from the \code{allD} data frame resulting from the database combination process. The original interactions data frame must have at least the following columns: the "K.ID" column (kinase Uniprot ID) and the "S.ID" column (substrate Uniprot ID), the "SID" column (unique character string identifiers for interactions), the "S.cc" column (substrateUPID.res.pos)
}
  \item{nK}{
a single string ("all","no", "drugs2data", "data") specifying whether the predicted edges should be included at all levels, none, only in the network that connects drug targets to data kinases or only between data sites and their kinases
}
}
\details{
If you add data to our network beyond the automatically provided one, the column "SID" needs to be filled, as any row with NA in that column will be removed (the weight correction depends on unique identifiers in that column).
}
\value{
This function returns a new \code{PKNlist} object that contains all interactions linking drug targets and drug targets to perturbed sites, and additional edges for integrator nodes (with SID="i…"), as well as all other network elements needed for optimisation
}
\author{
C. Terfve
}
