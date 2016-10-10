\name{BMsim}
\alias{BMsim}
\title{
Create and simulate a model
}
\description{
This function creates a pruned model and the simulation matrix associated with it
}
\usage{
BMsim(BM, nodes2link, BM.intAnd, pknList, targetsOn)
}
\arguments{
  \item{BM}{
a model data frame, i.e. a subset of a background model data frame (\code{interactions(pknList)}), that only contains those edges that are present in this particular model
}
  \item{BM.intAnd}{
a logical vector of the same length as the number of intermediates, that has the intermediates AND status associated with the model at hand 
}
  \item{nodes2link}{
a character vector of names of those nodes that we are trying to reach (i.e. they match nodes in the S.cc column of the background network, and should normally be data nodes, although they don't have to be)
}
  \item{pknList}{
A \code{PKNlist} object
}
  \item{targetsOn}{
a list of drugs targets, matching to the conditions in the data as stored in the \code{GMMbyCond} object (i.e. one element of this list contains the K/P affected under one condition, which itself can match to multiple drug treatments with identical targets).
}
}
\value{
a list with two elements
\item{BM.nopathR}{a model data frame with edges that cannot be connected to a data node removed}
\item{BM.sim}{a simulation matrix that has a row for each node and a column for each condition (set of perturbed kinases, possibly matching to multiple experimental treatments with the same affected K/P), x.ij=0 if node i is not perturbed under condition j, and 1 if it is perturbed. }

}
\author{
C. Terfve
}
