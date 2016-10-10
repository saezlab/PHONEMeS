\name{integratorWeights}
\alias{integratorWeights}
\title{
Integrator weights correction
}
\description{
Computes the correction to edges weight for integrator nodes, based on a set of solution models
}
\usage{
integratorWeights(models.sel=models.sel, pknList=pknList, models.sel.intgAnd=models.sel.intgAnd,models.sel.comb=models.sel.comb,optParam=optParam)
}
\arguments{
  \item{models.sel}{
list of solution networks
}
  \item{pknList}{
a \code{PKNlist} object that has at least elements: \code{integrators} (a character vector of integrator nodes names), \code{interactions} (the starting, complete, interactions data frame)
}
  \item{models.sel.intgAnd}{
matrix of combined intgAnd vectors for a set of solution vectors - this is only necessary if intgAsintm is TRUE, otherwise it can be set to NA
}
  \item{models.sel.comb}{
combined data frame of solution models - this is only necessary if \code{intgAsintm} is TRUE, otherwise it can be set to NA
}
  \item{optParam}{
an \code{OptParam} object with at least the elements \code{intgAsintm} (logical, should the integrators be treated as intermediates) and \code{cap} (numeric, maximum number of interactions that get copied)
}
}
\value{
a list with elements \code{integrators2paste} which is a data frame of interactions to be copied, and 2 vectors \code{AndBinF} and \code{NoneBinF} with an element for each integrator node. \code{AndBinF} is a new probability of the AND being sampled (if \code{intgAsintm}=TRUE) or an integer that corrects the bin size for the AND sampling (if \code{intgAsintm}=FALSE). \code{NoneBinF} is an integer that corrects the bin size for "no edge" (if \code{intgAsintm}=FALSE).
}
\author{
C. Terfve
}
