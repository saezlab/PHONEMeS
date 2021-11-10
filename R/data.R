#' The main PKN used in Phonemes to link deregulated phosphosites to perturbed kinases.
#' It combines kinase-substrate and protein-protein interactions.
#'
#'
#' @format The full PKN contains 56,738 interactions: 36,915 kinase-substrate interactions (1,688 kinases; 15,385 phosphosites); 19,823 protein-protein interactions
#' \describe{
#'     \item{source}{upstream protein}
#'     \item{interaction}{activation/inhibition (1/-1)}
#'     \item{target}{downstream protein/phosphorylation site}
#' }
#' @keywords PKN
#' @name phonemesPKN
#' @examples data("phonemesPKN")
#' @source \url{https://github.com/saezlab/PHONEMeS/blob/master/data-raw/phonemesPKN.R}
NULL

#' The KSN used as prior knowledge network when running vanilla Phonemes.
#' The KSN contains kinase-substrate interactions via a phosphorylation site
#'
#' @format The full PKN contains 37,631 kinase-substrate interactions via a phosphorylation site (1,710 kinases; 15,413 phosphosites)
#' \describe{
#'     \item{source}{upstream kinase/phosphorylation site}
#'     \item{interaction}{connection (1)}
#'     \item{target}{downstream phosphorylation site/kinase}
#' }
#' @keywords KSN
#' @name phonemesKSN
#' @examples data("phonemesKSN")
#' @source \url{https://github.com/saezlab/PHONEMeS/blob/master/data-raw/phonemesKSN.R}
NULL
