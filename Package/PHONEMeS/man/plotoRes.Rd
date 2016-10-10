\name{plotoRes}
\alias{plotoRes}
\title{
Produces optimisation diagnostic plots
}
\description{
This function takes a combined optimisation list as produced by oRes and plots a number of diagnostic curves to file
}
\usage{
plotoRes(opres, data.On=data.P, optParam=optParam)
}
\arguments{
  \item{opres}{
a list produced by \code{oRes}
}
  \item{data.On}{
a \code{GMMbyCond} object, as created by dataBycond (as it was fed into the optimisation)
}
  \item{optParam}{
an \code{OptParam} object with at least elements \code{resN} (the index of the optimisation, this is the suffix of the Resultsx.RData files, it will be used to name the pdf file "presN.pdf")
}
}
\details{
The plots that are produced are: 1. the average of the population and best model size, at each generation; 2.the drug specific scores of the best model of each generation, along generations; 3. the average drug specific scores of the population, along generations; 4. the total best model and population average scores along generations (without size penalty); 5. the total best model and population average scores along generations (with size penalty); 6.the evolution of frequency of edges along generations (1 curve for each edge), separated in 3 categories 6.A. edges that have a final frequency between 0.4 and 0.6, 6.B. edges that have a final frequency above 0.6, and 6.C. edges that have a final frequency below 0.4.  Note that all these numbers are computed based on a population size of 5000.
}
\value{
This function doesn't return anything, it just produces plots that are stored in a single pdf file.
}

\author{
C. Terfve
}
