#' Run PHONEMeS
#'
#' This function runs CARNIVAL with the input of phosphoproteomic data (phosphosites and kinases).
#' The prior knowledge network used is the combination of protein-protein and protein-phosphosite
#' interactions from omnipath. Before running CARNIVAL the network is pruned by removing nodes n_steps
#' upstream and downstream of measurements and inputs, respectively.
#'
#' @param inputObj named vector of perturbation targets. Either 1 (up regulated) or -1 (down regulated)
#' @param measObj named vector of the measurements
#' @param netObj data frame of the prior knowledge network
#' @param rmNodes character vector of nodes to remove from prior knowledge network
#' @param pruning logic, set to TRUE if network should be pruned (recommended)
#' @param n_steps_pruning integer giving the order of the neighborhood
#' @param options An object of type \dQuote{\code{list}} defining the run parameters CARNIVAL in PHONEMeS.
#'   Use the \code{\link{default_phonemes_options}} function to create a list with default parameter settings.
#'   If cplex or cbc are chosen as the solver, the parameter solverPath needs to be supplied.
#' @return List of CARNIVAL results and final inputObj, measObj, netObj used
#' @importFrom dplyr %>%
#' @export

run_phonemes <- function(inputObj,
                         measObj,
                         netObj = phonemesPKN,
                         rmNodes = NULL,
                         pruning = TRUE,
                         n_steps_pruning = 50,
                         carnival_options){

  netObj <- netObj %>% dplyr::filter(!(source %in% rmNodes | target %in% rmNodes))

  # Remove input and measurements not part of the PKN
  inputObj <- inputObj[names(inputObj) %in% netObj$source]
  measObj <- measObj[names(measObj) %in% netObj$target]

  if (pruning) {
    # Remove nodes n_steps downstream of perturbations
    meta_g <- igraph::graph_from_data_frame(netObj[,c("source","target",'interaction')],directed = TRUE)
    dn_nbours <- igraph::ego(graph = meta_g, order = n_steps_pruning, nodes = names(inputObj), mode = "out")
    sub_nodes <- c(unique(names(unlist(dn_nbours))), names(inputObj))
    pruned_PKN <- netObj %>% dplyr::filter(source %in% sub_nodes & target %in% sub_nodes)

    # Remove nodes n_steps upstream of perturbations
    up_nbours <- igraph::ego(graph = meta_g, order = n_steps_pruning, nodes = names(measObj), mode = "in")
    up_nodes <- c(unique(names(unlist(up_nbours))), names(measObj))
    netObj <- pruned_PKN %>% dplyr::filter(source %in% up_nodes & target %in% up_nodes)
  }

  # Remove input and measurements not part of the PKN (2nd pruning because some nodes have disapeared)
  inputObj <- inputObj[names(inputObj) %in% netObj$source]
  measObj <- measObj[names(measObj) %in% netObj$target]

  message(paste("Input nodes:", length(inputObj),
                "\nMeasurement nodes:", length(measObj),
                "\nNetwork nodes:", length(unique(c(netObj$source, netObj$target))),
                "\nNetwork edges:", nrow(netObj)))

  check_carnival_options(carnival_options)
  resCarnival <- CARNIVAL::runVanillaCarnival(perturbations = inputObj,
                                              measurements = measObj,
                                              priorKnowledgeNetwork = netObj,
                                              carnivalOptions = carnival_options)

  # Remove nodes with 0 weight
  resCarnival$weightedSIF <- resCarnival$weightedSIF %>% dplyr::filter(Weight != 0)
  resCarnival$nodesAttributes <- resCarnival$nodesAttributes %>% dplyr::filter(Node %in% union(resCarnival$weightedSIF$Node1, resCarnival$weightedSIF$Node2))

  # Add degree to attributes
  in_degree <- resCarnival$weightedSIF %>% dplyr::group_by(Node2) %>%
    dplyr::summarise(in_degree  = dplyr::n()) %>%
    dplyr::rename(Node = "Node2")

  out_degree <- resCarnival$weightedSIF %>%
    dplyr::group_by(Node1) %>%
    dplyr::summarise(out_degree  = dplyr::n()) %>%
    dplyr::rename(Node = "Node1")

  degree_df <- base::merge(in_degree, out_degree, by = "Node", all = TRUE) %>%
    as.data.frame() %>%
    tidyr::replace_na(list(in_degree = 0, out_degree = 0)) %>%
    dplyr::mutate(tot_degree = rowSums(dplyr::across(c(in_degree, out_degree))))

  resCarnival$nodesAttributes <- base::merge(resCarnival$nodesAttributes,degree_df, by = "Node", all = TRUE) %>% as.data.frame()

  return(list(res = resCarnival,
              network = netObj,
              measurements = measObj,
              inputs = inputObj))

}

#' Reattach_psites
#'
#' This function readd links between phosphosite and their correpsonding proteins
#'
#' @param phonemes_res phonemes result from the run_phonemes function
#' @return List of PHONEMES results and final inputObj, measObj, netObj used, with psites attached
#' @export
#'
reattach_psites <- function(phonemes_res)
{
  sif <- phonemes_res$res$weightedSIF
  att <- phonemes_res$res$nodesAttributes

  phospho_prots <- data.frame(sif[grepl("_",sif$Node2),3])
  names(phospho_prots) <- "Node1"
  phospho_prots$Node2 <- gsub("_.*","",phospho_prots$Node1)
  phospho_prots$Sign <- 1
  phospho_prots$Weight <- 1
  phospho_prots <- phospho_prots[phospho_prots$Node2 %in% att$Node,]

  if(length(phospho_prots[,1]) > 0)
  {
    sif <- as.data.frame(rbind(sif, phospho_prots))
    sif <- unique(sif)
  } else
  {
    print("No psites to attach")
  }

  phonemes_res$res$weightedSIF <- sif
  phonemes_res$res$nodesAttributes <- att

  return(phonemes_res)
}

#' get_protein_network
#'
#' This function readd links between phosphosite and their correpsonding proteins
#'
#' @param phonemes_res phonemes result from the run_phonemes function
#' @return Phonemes network only consisting of protein-protein interactions
#' @importFrom dplyr %>%
#' @export
#'
get_protein_network <- function(phonemes_res)
{
  sif <- phonemes_res$res$weightedSIF
  att <- phonemes_res$res$nodesAttributes

  sif <- sif %>% dplyr::filter(!grepl(pattern = "[a-zA-Z0-9]_[a-zA-Z0-9]", sif$Node2))
  att <- att %>% dplyr::filter(Node %in% union(sif$Node1, sif$Node2))

  # Add protein degree to attributes
  in_degree <- sif %>% dplyr::group_by(Node2) %>%
    dplyr::summarise(protein_in_degree  = dplyr::n()) %>%
    dplyr::rename(Node = "Node2")

  out_degree <- sif %>%
    dplyr::group_by(Node1) %>%
    dplyr::summarise(protein_out_degree  = dplyr::n()) %>%
    dplyr::rename(Node = "Node1")

  degree_df <- base::merge(in_degree, out_degree, by = "Node", all = TRUE) %>%
    as.data.frame() %>%
    tidyr::replace_na(list(protein_in_degree = 0, protein_out_degree = 0)) %>%
    dplyr::mutate(protein_tot_degree = rowSums(dplyr::across(c(protein_in_degree, protein_out_degree))))

  att <- base::merge(att,degree_df, by = "Node", all = TRUE) %>% as.data.frame()

  protein_network <- list(weightedSIF = sif,
                          nodesAttributes = att)

  return(protein_network)
}

#' extract_subnetwork
#'
#' This function extracts smaller sub networks from the run_phonemes output
#'
#' @param phonemes_res Phonemes result from the run_phonemes function
#' @param targets Network nodes, starting point for the extraction of the sub network
#' @param n_steps Number of steps to extract down- or upstream of targets
#' @param mode Character constant to specify direction of the extraction. "In" for upstream nodes, "out" for downstream nodes and "all" for both.
#' @return Phonemes sub network
#' @importFrom dplyr %>%
#' @export
#'
extract_subnetwork <- function(phonemes_res, targets, n_steps = 3, mode = "all")
{
  sif <- phonemes_res$res$weightedSIF
  att <- phonemes_res$res$nodesAttributes

  targets <- targets[targets %in% att$Node]
  if(purrr::is_empty(targets)){
    warning("No target found in network. Returning empty data.frame")
    return(list(weightedSIF = data.frame(),
                nodesAttributes = data.frame()))
  }

  meta_g <- igraph::graph_from_data_frame(sif[,c("Node1","Node2",'Sign')],directed = TRUE)

  if (mode %in% c("in", "out")) {
    dn_nbours <- igraph::ego(graph = meta_g, order = n_steps, nodes = targets, mode = mode)
    sub_nodes <- c(unique(names(unlist(dn_nbours))), targets)
  } else if (mode %in% "all") {
    dn_nbours_in <- igraph::ego(graph = meta_g, order = n_steps, nodes = targets, mode = "in")
    dn_nbours_out <- igraph::ego(graph = meta_g, order = n_steps, nodes = targets, mode = "out")
    sub_nodes <- c(unique(c(names(unlist(dn_nbours_in)), names(unlist(dn_nbours_out)))), targets)
  }

  sif <- sif %>% dplyr::filter(Node1 %in% sub_nodes & Node2 %in% sub_nodes)
  att <- att %>% dplyr::filter(Node %in% union(sif$Node1, sif$Node2))

  subnetwork <- list(weightedSIF = sif,
                     nodesAttributes = att)

  return(subnetwork)

}
