#' The PKN used in Phonemes to link deregulated phosphosites to perturbed kinases.
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
