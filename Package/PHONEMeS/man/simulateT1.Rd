\name{simulateT1}
\alias{simulateT1}
\title{
One update step in simulation
}
\description{
This function does one update step in a simulation process
}
\usage{
simulateT1(nodes.t=nodes.t,intermediates.AND=intermediates.AND,integrators.AND=integrators.AND,intm.inN=intm.inN,intg.inN=intg.inN,sink.inN=sink.inN)
}
\arguments{
  \item{nodes.t}{
numeric matrix with two columns being the state of each node at state t and state t+1, status of nodes are 0 or 1
}
  \item{intermediates.AND}{
for each intermediate in the network, this is a logical that says if inputs enter it with an AND gate or not
}
  \item{integrators.AND}{
for each integrator in the network, this is a logical that says if inputs enter it with an AND gate or not. This is NA if the integrators are not considered as intermediates
}
  \item{intm.inN}{
a list with an entry for each intermediate that has inputs, containing a character vector with the names of input nodes of the node
}
  \item{intg.inN}{
a list with an entry for each integrator that has inputs, containing a character vector with the names of input nodes of the node
}
  \item{sink.inN}{
a list with an entry for each sink that has inputs, containing a character vector with the names of input nodes of the node
}
}
\details{
This function computes the status of each node (recorded in column 2 of nodes.T) at t+1 as a function of the status of its parents at t (column 1) according to a model
}
\value{
the updated nodes.t matrix
}
\author{
C. Terfve
}
