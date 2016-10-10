\name{PKNlist}
\alias{PKNlist}
\title{
Create a \code{PKNlist} object
}
\description{
This function creates a \code{PKNlist} object which contains all the elements based on a background network as created by \code{buildPKN}, that are needed for the optimisation
}
\usage{
PKNlist(PKN, targets.On, data.On)
}
\arguments{
  \item{PKN}{
a \code{KPSbg} object created by \code{buildPKN}
}
  \item{targets.On}{
a list of drug targets (as character vectors) by condition, matching to the \code{data.On} object
}
  \item{data.On}{
a \code{GMMbyCond} object as created by \code{dataBycond}
}
}
\value{
a \code{PKNlist} object with elements \code{interactions} (a condensed form of the version in the original network in the object created by \code{buildPKN}, that aggregates identical interactions from multiple sources, with weights modified accordingly), \code{interactionsD} (a subset of interactions that only contains experimentally supported direct interactions between data points and drug targets), \code{intermediates} (a character vector of intermediate nodes), \code{sinks} (a character vector of sink nodes), \code{integrators} (a character vector of integrator nodes).
}
\author{
C. Terfve
}
