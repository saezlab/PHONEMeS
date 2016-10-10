\name{pathsD}
\alias{pathsD}
\title{
Max scoring path
}
\description{
Extracts, from a background network with optimised weights, the maximum scoring paths from a kinase to a set of data nodes
}
\usage{
pathsD(nwTable, nM,nodes2link,drug)
}
\arguments{
  \item{nwTable}{
a background network data frame (from \code{interactions(pknList)}), with appropriate optimised weights in the ntag column (this is what is used to decide on the max input)
}
  \item{nM}{
the number of models on which the frequency is based}
  \item{nodes2link}{
a character vector of data nodes that we want to connect
}
  \item{drug}{
a single character string, containing the name of the K/P that we want to connect (a single UPID)
}
}
\details{
The maximum input edge is picked based on the "ntag" column of the \code{nwTable} data frame, so it should contain \code{opres$FE[match(nwTable$SID, rownames(opres$FE)),dim(opres$FE)[2]]} where \code{opres} is the list resulting from the \code{oRes} function. If we are working with combined optimisations, then the \code{ntag} column should contained the averaged equivalents of the vector just mentioned. 
}
\value{
a data frame of interactions that is a subset of the starting one
}

\author{
C. Terfve
}
