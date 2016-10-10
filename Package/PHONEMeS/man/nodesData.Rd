\name{nodesData}
\alias{nodesData}
\title{
Nodes On/Off in a network, given some data
}
\description{
This function extracts the nodes that are On/Off in a network, according to a given dataset
}
\usage{
nodesData(data.On, dataGMM, pknList)
}
\arguments{
  \item{data.On}{
a \code{GMMbyCond} object as created by \code{dataBycond}
}
  \item{dataGMM}{
a \code{GMMres} object
}
  \item{pknList}{
a \code{PKNlist} object
}
}

\value{
This function returns a list that contains elements \code{On}, \code{Off}, \code{Onlist} and \code{Offlist}, which contain character vectors of names of nodes that are On/Off in the network given by \code{pknList}, according to the data in \code{data.On} and \code{dataGMM}, either as single character vectors (\code{On}, \code{Off}) or organised by condition and treatments like in \code{data.On} (\code{Onlist}, \code{Offlist}).
}
\author{
C. Terfve
}
