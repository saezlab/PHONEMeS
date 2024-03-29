% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vanilla_phonemes.R
\name{run_vanilla_phonemes}
\alias{run_vanilla_phonemes}
\title{Run vanilla PHONEMeS (PHONEMeS v1.0.0)}
\usage{
run_vanilla_phonemes(
  inputObj,
  measObj,
  netObj = phonemesKSN,
  rmNodes = NULL,
  pruning = FALSE,
  n_steps_pruning = 50,
  carnival_options
)
}
\arguments{
\item{inputObj}{named vector of perturbation targets}

\item{measObj}{named vector of the measurements}

\item{netObj}{data frame of the prior knowledge network}

\item{rmNodes}{character vector of nodes to remove from prior knowledge network}

\item{pruning}{logic, set to TRUE if network should be pruned (default = FALSE)}

\item{n_steps_pruning}{integer giving the order of the neighborhood}

\item{options}{An object of type \dQuote{\code{list}} defining the run parameters CARNIVAL in PHONEMeS.
Use the \code{\link{default_phonemes_options}} function to create a list with default parameter settings.
If cplex or cbc are chosen as the solver, the parameter solverPath needs to be supplied.}
}
\value{
List of CARNIVAL results and final inputObj, measObj, netObj used
}
\description{
The code for this function is comparable to the orginal PHONEMeS package. An unsigned
kinase-substrate network is used as prior knowledge connecting perturbed kinases with phosphosites.
}
