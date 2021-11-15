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
#' @param solverPath path to the solver
#' @param solver one of the solvers available from getSupportedSolvers()
#' @param timelimit solver time limit in seconds
#' @param mipGAP CPLEX parameter: absolute tolerance on the gap
#' @param poolrelGAP CPLEX/Cbc parameter: Allowed relative gap of accepted
#' @return List of CARNIVAL results and final inputObj, measObj, netObj used
#' @importFrom dplyr %>%
#' @export

run_phonemes <- function(inputObj,
                          measObj,
                          netObj = phonemesPKN,
                          rmNodes = NULL,
                          pruning = TRUE,
                          n_steps_pruning = 50,
                          solverPath,
                          solver = "cplex",
                          timelimit = 7200,
                          mipGAP = 0.05,
                          poolrelGAP = 0){

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

  # Turn inputObj and measObj into data.frames for Carnival 1.3.0
  inputObj <-  as.data.frame(t(inputObj))
  measObj <-  as.data.frame(t(measObj))

  message(paste("Input nodes:", length(inputObj),
                "\nMeasurement nodes:", length(measObj),
                "\nNetwork nodes:", length(unique(c(netObj$source, netObj$target))),
                "\nNetwork edges:", nrow(netObj)))

  resCarnival <- CARNIVAL::runCARNIVAL(inputObj = inputObj,
                                      measObj = measObj,
                                      netObj = netObj,
                                      solverPath = solverPath,
                                      solver = solver,
                                      timelimit = timelimit,
                                      mipGAP = mipGAP,
                                      poolrelGAP = poolrelGAP)

  # Remove nodes from Carnival results that have no up- or downAct
  resCarnival$nodesAttributes <- as.data.frame(resCarnival$nodesAttributes) %>% dplyr::mutate(dplyr::across(c(ZeroAct, UpAct, DownAct, AvgAct), as.double))
  zeroNodes <- resCarnival$nodesAttributes %>% dplyr::filter(UpAct == 0 & DownAct == 0) %>% dplyr::pull(Node)

  resCarnival$nodesAttributes <- resCarnival$nodesAttributes %>% dplyr::filter(!Node %in% zeroNodes)
  rm(zeroNodes)

  resCarnival$weightedSIF <- as.data.frame(resCarnival$weightedSIF) %>%
    dplyr::mutate(dplyr::across(Weight, as.double)) %>%
    dplyr::filter(Node1 %in% resCarnival$nodesAttributes$Node & Node2 %in% resCarnival$nodesAttributes$Node)

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
