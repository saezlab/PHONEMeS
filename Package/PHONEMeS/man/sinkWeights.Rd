\name{sinkWeights}
\alias{sinkWeights}
\title{
Sinks weight correction
}
\description{
Computes the correction to edges weight for sink nodes, based on a set of solution models
}
\usage{
sinkWeights(pknList=pknList, optParam=optParam, models.sel.comb=models.sel.comb, models.sel=models.sel)
}
\arguments{
  \item{optParam}{
an \code{OptParam} object with at least the element \code{cap} (numeric, maximum number of interactions that get copied)
}
  \item{pknList}{
a \code{PKNlist} object that has at least elements: \code{sinks} (a character vector of sinks nodes names), \code{interactions} (the starting, complete, interactions data frame) 
}
  \item{models.sel.comb}{
data frame of combined solution networks
}
  \item{models.sel}{
list of solution networks
}
}
\value{
a data frame of interactions to be copied
}
\author{
C. Terfve
}