\name{run1model}
\alias{run1model}
\title{
Creates and runs 1 model
}
\description{
This function creates a model, simulates it and outputs its score, given a data set
}
\usage{
run1model(pknList=pknList, targetsOn=targets.P,dataOn=data.P,FlipP=rep(0.5, length(intermediates(pknList))),AndBinF=rep(0, length(integrators(pknList))),NoneBinF=rep(0, length(integrators(pknList))),dataGMM=dataGMM,optParam)
}
\arguments{
  \item{pknList}{
a \code{PKNlist} object
}
  \item{optParam}{
\code{cstart} (logical, should the \code{interactionsD} set of interactions be added to every model, normally only set to TRUE at generation 1), and \code{intgAsintm} (logical, should integrators be treated like intermediates)
}
  \item{targetsOn}{
a list of string vectors containing the UPIDs of kinases hit by each set of drugs
}
  \item{dataOn}{
a \code{GMMbyCond} object
}
  \item{FlipP}{
a numeric vector, the probability of AND for each intermediate
}
  \item{AndBinF}{
a numeric vector, either an integer giving the number of bins for the AND gate of each integrator (if \code{intgAsintm(optParam)=FALSE}), or the probability of AND for each integrator (if \code{intgAsintm(optParam)=TRUE})
}
  \item{NoneBinF}{
a numeric vector (of the same length as the number of integrator nodes), the number of bins for the no-gate option, for each integrator node
}
  \item{dataGMM}{
a GMMres object
}
}
\value{
a list with elements: model (the model, as created by \code{createModel}, \code{score} (the score of the model), \code{scoresList} (a list of scores in the same format as \code{dataOn} but with a single numeric value in each element, which is the score produced by the data under that drug), \code{intAnd/intgAnd} (a logical vector for each intermediate/integrator: AND or OR gated). \code{AndBinF} and \code{NoneBinF} are NA if \code{intgAsintm} is TRUE (deprecated).
}
\author{
C. Terfve
}
