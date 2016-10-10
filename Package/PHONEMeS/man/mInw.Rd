\name{mInw}
\alias{mInw}
\title{
Create maximum input network
}
\description{
Extracts, from a background network with optimised weights, the maximum input network with tolerance
}
\usage{
mInw(nwTable, intgNone, nodes2link, tol, targets.On)
}
\arguments{
  \item{nwTable}{
a background network data frame (from \code{interactions(pknList)}), with appropriate optimised weights in the ntag column (this is what is used to decide on the max input)
}
  \item{intgNone}{
a numeric vector of the same length as the number of integrators, that has the "no input" weight for each integrator; this vector should be named (mapping to the integrators names)
}
  \item{nodesOF}{
a list created by \code{nodesData}, containing an element "On" which is a character vector of names of those nodes that we are trying to reach (i.e. they match nodes in the S.cc column of the background network, and should normally be data nodes, although they don't have to be), and an element "OnList" containing the same but organised as a list by condition
}
  \item{tol}{
a single number between 0 and 1, the relative tolerance around the maximum input edge that we allow (!not to be confused with the optimisation tolerance that is the tolerance used in selecting best models for weights correction at each generation)
}
  \item{targets.On}{
a list of character vectors containing the K/P perturbed under each condition (one condition=one element of the list)
}
}
\details{
The maximum input edge is picked based on the "ntag" column of the \code{nwTable} data frame, so it should contain the appropriate values: if we want a standard ntag network, then ntag should be \code{interactions(pknList)} (of the last generation), if we want a frequency based network, then \code{nwTable$ntag} should be \code{opres$FE[match(nwTable$SID, rownames(opres$FE)),dim(opres$FE)[2]]} where \code{opres} is the list resulting from the \code{oRes} function. If we are working with combined optimisations, then the \code{ntag} column should contained the averaged equivalents of the vectors just mentioned. 
}
\value{
a data frame of interactions that is a subset of the starting one
}

\author{
C. Terfve
}
