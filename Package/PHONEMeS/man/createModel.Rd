\name{createModel}
\alias{createModel}
\title{
Sample a candidate model
}
\description{
This function samples an interaction object (\code{PKNlist}) to make a candidate model
}
\usage{
createModel(pknList=pknList,AndBinF=rep(1, length(integrators(pknList))),NoneBinF=rep(1, length(integrators(pknList))),optParam=optParam)
}
\arguments{
  \item{pknList}{
a \code{PKNlist} object produced by \code{buildPKN}, with elements \code{interactions} (an interaction data frame), \code{sinks} (a character vector of sink nodes), \code{species} (a character vector of nodes), \code{interactionsD} (a subset of interactions that only contains experimentally supported direct interactions between drug targets and data points), \code{intermediates} (a character vector of intermediate nodes), \code{integrators} (a character vector of integrator nodes)
}
  \item{AndBinF}{
a numeric vector of the same length as the number of integrators
}
  \item{NoneBinF}{
a numeric vector of the same length as the number of integrators
}
  \item{optParam}{
an \code{OptParam} object with at least the elements \code{intgAsintm} (logical, TRUE if we want integrators to be sampled the same way intermediates are) and \code{cstart} (logical, should the interactionsD set of interactions be added to every model)
}
}
\value{
a data frame of interactions that is a subset of the interactions in \code{interactions(pknList)}
}
\author{
C. Terfve
}
