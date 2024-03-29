% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_phonemes.R
\name{extract_subnetwork}
\alias{extract_subnetwork}
\title{extract_subnetwork}
\usage{
extract_subnetwork(phonemes_res, targets, n_steps = 3, mode = "all")
}
\arguments{
\item{phonemes_res}{Phonemes result from the run_phonemes function}

\item{targets}{Network nodes, starting point for the extraction of the sub network}

\item{n_steps}{Number of steps to extract down- or upstream of targets}

\item{mode}{Character constant to specify direction of the extraction. "In" for upstream nodes, "out" for downstream nodes and "all" for both.}
}
\value{
Phonemes sub network
}
\description{
This function extracts smaller sub networks from the run_phonemes output
}
