#' regulatory_psites
#'
#' This function  extract regulatory phosphites from the phonemes network, e.i.
#' psites that are differentially regulated on proteins that are found in
#' the phonemes network
#'
#' @param phonemes_res carnival result from the run_carnival function
#' @param th_act threshold for node activity to include in the benchmark
#' @return list with two elements: the regulatory psites with predicted mode of
#' regulations and the kinases that catalyse their phosphorilation
#' @export
#'
regulatory_psites <- function(phonemes_res, th_act = 0)
{
  sif <- as.data.frame(phonemes_res$res$weightedSIF)
  att <- as.data.frame(phonemes_res$res$nodesAttributes)
  att <- att[att$AvgAct != 0,]

  predicted_phospho_function <- sif[grepl("_",sif$Node1),]
  predicted_phospho_function <- sif[grepl("_",sif$Node1) | sif$Node2 %in% predicted_phospho_function$Node1,]
  predicted_phospho_function <- predicted_phospho_function[
    grepl("_",predicted_phospho_function$Node1) |
      predicted_phospho_function$Node1 != gsub("_.*","",predicted_phospho_function$Node2),
  ]
  predicted_phospho_function <- predicted_phospho_function[
    !grepl("_",predicted_phospho_function$Node1) |
      predicted_phospho_function$Node1 %in% predicted_phospho_function$Node2,
  ]
  functionnal_kinases <- predicted_phospho_function[
    !grepl("_",predicted_phospho_function$Node1),]

  predicted_phospho_function <- predicted_phospho_function[
    grepl("_",predicted_phospho_function$Node1),]

  predicted_phospho_function$Sign <- sapply(predicted_phospho_function$Node1, function(x, att){
    protein <- gsub("_.*","",x)
    if(abs(as.numeric(att[att$Node == protein,"AvgAct"])) >= th_act)
    {
      mor <- sign(as.numeric(att[att$Node == protein,"AvgAct"]) * as.numeric(att[att$Node == x, "AvgAct"]))
    } else
    {
      mor <- NA
    }

    return(mor)
  }, att = att)

  predicted_phospho_function <- predicted_phospho_function[,c(1,2)]
  names(predicted_phospho_function) <- c("psite","mor")

  return(list("regulatory_psites" = predicted_phospho_function, "functional_kinases" = functionnal_kinases))
}

#' benchmark_regulatory_psites
#'
#' This function check the coherence of regulatory phosphosites predicted by
#' phonemes with ground truth from phosphositeplus
#'
#' @param sif regulatory_psites, e.i. the first element of the list returned by
#' regulatory_psites function
#' @return a table comparing the predicted regulatory psites with the regulatory
#' psites of phosphositeplus
#' @export
#'
benchmark_regulatory_psites <- function(regulatory_psites)
{
  regulatory_psites <- merge(regulatory_psites, reg_sites, by = "psite",all.x = T)

  regulatory_psites_complete <- regulatory_psites[
    complete.cases(regulatory_psites),]

  return(regulatory_psites_complete)
}

#' process_phosphositeplus_regpsites
#'
#' This function process the tab file of regulatory psites of phosphositeplus,
#' downloaded from https://www.phosphosite.org/staticDownloads
#'
#' @param regpsite_file path to the tab file of regulatory psites of phosphositeplus
#' @return a formatted table of regulatory psites of phosphositeplus
#' with their known mode of regulation
#' @export
#'
process_phosphositeplus_regpsites <- function(ppsp_regpsites_file_path)
{
  reg_sites_full <- as.data.frame(
    read_delim(ppsp_regpsites_file_path, "\t", escape_double = FALSE, trim_ws = TRUE))
  reg_sites_full <- reg_sites_full[grepl("human",reg_sites_full$ORGANISM),]
  reg_sites_full$psite <- paste(reg_sites_full$GENE, gsub("-.*","",reg_sites_full$MOD_RSD), sep = "_")
  reg_sites <- reg_sites_full[,c(22,12)]
  reg_sites$ppsp_mor <- ifelse(grepl("induced",reg_sites$ON_FUNCTION),
                               1,
                               ifelse(grepl("inhib",reg_sites$ON_FUNCTION),-1,NA))
  reg_sites <- reg_sites[complete.cases(reg_sites),]
  reg_sites <- reg_sites[,-2]
  return(reg_sites)
}

