#' Run vanilla PHONEMeS (PHONEMeS v1.0.0)
#'
#' The code for this function is comparable to the orginal PHONEMeS package. An unsigned
#' kinase-substrate network is used as prior knowledge connecting perturbed kinases with phosphosites.
#'
#' @param inputObj named vector of perturbation targets. Omit to run InvCARNIVAL
#' @param measObj named vector of the measurements
#' @param netObj data frame of the prior knowledge network
#' @param rmNodes character vector of nodes to remove from prior knowledge network
#' @param pruning logic, set to TRUE if network should be pruned (default = FALSE)
#' @param n_steps_pruning integer giving the order of the neighborhood
#' @param solverPath path to the solver
#' @param solver one of the solvers available from getSupportedSolvers()
#' @param timelimit solver time limit in seconds
#' @param mipGAP CPLEX parameter: absolute tolerance on the gap
#' @param poolrelGAP CPLEX/Cbc parameter: Allowed relative gap of accepted
#' @return List of CARNIVAL results and final inputObj, measObj, netObj used
#' @importFrom dplyr %>%
#' @export

run_vanilla_phonemes <- function(inputObj = NULL,
                                 measObj,
                                 netObj = phonemesKSN,
                                 rmNodes = NULL,
                                 pruning = FALSE,
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
    meta_g <- igraph::graph_from_data_frame(netObj[, c("source", "target", 'interaction')], directed = TRUE)

    if (!is.null(inputObj)) {
      # Remove nodes n_steps downstream of perturbations
      dn_nbours <- igraph::ego(graph = meta_g, order = n_steps_pruning, nodes = names(inputObj), mode = "out")
      sub_nodes <- c(unique(names(unlist(dn_nbours))), names(inputObj))
      netObj <- netObj %>% dplyr::filter(source %in% sub_nodes & target %in% sub_nodes)
    }
    # Remove nodes n_steps upstream of perturbations
    up_nbours <- igraph::ego(graph = meta_g, order = n_steps_pruning, nodes = names(measObj), mode = "in")
    up_nodes <- c(unique(names(unlist(up_nbours))), names(measObj))
    netObj <- netObj %>% dplyr::filter(source %in% up_nodes & target %in% up_nodes)
  }

  # Remove input and measurements not part of the PKN (2nd pruning because some nodes have disapeared)
  inputObj <- inputObj[names(inputObj) %in% netObj$source] %>% abs()
  measObj <- measObj[names(measObj) %in% netObj$target] %>% abs()

  # Turn inputObj and measObj into data.frames for Carnival 1.3.0
  if (!is.null(inputObj)) {
    inputObj <- as.data.frame(t(inputObj))
  }
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

  resCarnival$nodesAttributes <- resCarnival$nodesAttributes %>%
    dplyr::filter(!Node %in% zeroNodes) %>%
    dplyr::select(-c(ZeroAct, UpAct, DownAct))

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
