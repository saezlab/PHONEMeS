\name{intermediateWeights}
\alias{intermediateWeights}
\title{
Correction to intermediate weights
}
\description{
Computes the correction to current intermediate weights given a set of solution networks
}
\usage{
intermediateWeights(pknList=pknList,models.sel.intAnd=models.sel.intAnd,models.sel=models.sel,models.sel.comb=models.sel.comb,optParam=optParam)
}
\arguments{
  \item{pknList}{
a \code{PKNlist} object that has at least elements: \code{intermediates} (a character vector of intermediate nodes names), \code{interactions} (the starting, complete, interactions data frame)
}
  \item{models.sel.intAnd}{
matrix of combined intAnd vectors for a set of solution vectors
}
  \item{models.sel}{
list of solution models
}
  \item{models.sel.comb}{
combined data frame of solution models
}
  \item{optParam}{
an \code{OptParam} object with at least the element \code{cap} (numeric, maximum number of interactions that get copied)
}
}
\value{
a list with two elements: a data frame of interactions that should be duplicated and a vector of new AndFlipF (the new probability, for each intermediate node, of the AND being sampled)
}
\author{
C. Terfve
}
