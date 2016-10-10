\name{oRes}
\alias{oRes}
\title{
Import and combine optimisation results
}
\description{
Function that imports and combines the results of an optimisation.
}
\usage{
oRes(data.On, optParam, pknList)
}
\arguments{
  \item{optParam}{
an \code{OptParam} object with at least elements \code{resN} (the index of the optimisation, this is the suffix of the Resultsx.RData files, and matches the suffix of the results folder \code{wd/Results_resN} where all results will be written, \code{tol} (the tolerance that was used for the optimisation, a number between 0 and 1 if relative tolerance, or a number larger or equal to 1 if absolute tolerance), \code{nG} (the number of generations that the optimisation was run for)
}
  \item{data.On}{
a \code{GMMbyCond} object, as created by dataBycond (as it was fed into the optimisation)
}
  \item{pknList}{
a \code{PKNlist} object as it was fed in the optimisation process
}
}
\details{
This function imports all the output and results files from the optimisation process and returns summarised objects that contain the evolution of scores, weights, etc.
}
\value{
This function returns a list of list: 

\item{nM}{vector, the number of models that pass the tolerance threshold, at each generation}
\item{sM}{matrix, score of the best model of the generation, by drug (rows=generations, cols=drugs, in the same order as in data.On)}
\item{sM.avg}{matrix, average score of the generation, by drug (rows=generations, cols=drugs, in the same order as in data.On)}
\item{sAll}{matrix, total score of each model, at each generation (rows=generations, cols=models in population, fixed to 5000 currently)}
\item{Msize}{matrix, size of all models, at each generation, mapping to sAll}
\item{BMGsize}{vector, size of the best model of each generation}
\item{BM.s}{vector, total score of the best model of the generation (with size penalty if any)}
\item{FE}{matrix, counts of each edge at each generation (rows=edges (in the same order as in c.I.list$complete.I), cols=generations); to get the frequency from these numbers, each entry should be divided by the number of models in the population, currently fixed to 5000}
\item{intgNone}{matrix, counts of "no inputs" for each integrator at every generation (rows=integrator nodes (in the same order as in c.I.list$integrators), cols=generations); to get the frequency from these numbers, each entry should be divided by the number of models in the population, currently fixed to 5000}
\item{intgAnd}{matrix, counts of "ands" for each integrator at every generation (rows=integrator nodes (in the same order as in c.I.list$integrators), cols=generations); to get the frequency from these numbers, each entry should be divided by the number of models in the population, currently fixed to 5000}
\item{G1.freq}{matrix, ntag of each edge at each generation (rows=edges (in the same order as in c.I.list$complete.I), cols=generations); these are the number of bins for each edge used for sampling at the next generation}
\item{G1.flipP}{matrix, probability of "and" for each intermediate at every generation (rows=intermediate nodes (in the same order as in c.I.list$intermediates), cols=generations)}
\item{G1.AndBinF}{matrix, ntag of the "AND" integrator hyperedge at each generation (rows=integrator nodes (in the same order as in c.I.list$integrators), cols=generations); these are the number of bins for the AND hyperedge for each integrator node,  used for sampling at the next generation}
\item{G1.NoneBinF}{matrix, ntag of the "none" integrator hyperedge at each generation (rows=integrator nodes (in the same order as in c.I.list$integrators), cols=generations); these are the number of bins for the "no input" edge for each integrator node,  used for sampling at the next generation}
}
\author{
C. Terfve
}