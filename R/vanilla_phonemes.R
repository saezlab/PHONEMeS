#' Run vanilla PHONEMeS (PHONEMeS v1.0.0)
#'
#' The code for this function is comparable to the orginal PHONEMeS package. An unsigned
#' kinase-substrate network is used as prior knowledge connecting perturbed kinases with phosphosites.
#'
#' @param inputObj named vector of perturbation targets
#' @param measObj named vector of the measurements
#' @param netObj data frame of the prior knowledge network
#' @param rmNodes character vector of nodes to remove from prior knowledge network
#' @param pruning logic, set to TRUE if network should be pruned (default = FALSE)
#' @param n_steps_pruning integer giving the order of the neighborhood
#' @param options An object of type \dQuote{\code{list}} defining the run parameters CARNIVAL in PHONEMeS.
#'   Use the \code{\link{default_phonemes_options}} function to create a list with default parameter settings.
#'   If cplex or cbc are chosen as the solver, the parameter solverPath needs to be supplied.
#' @return List of CARNIVAL results and final inputObj, measObj, netObj used
#' @importFrom dplyr %>%
#' @export

run_vanilla_phonemes <- function(inputObj,
                                 measObj,
                                 netObj = phonemesKSN,
                                 rmNodes = NULL,
                                 pruning = FALSE,
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
  inputObj <- inputObj[names(inputObj) %in% netObj$source] %>% abs()
  measObj <- measObj[names(measObj) %in% netObj$target] %>% abs()

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
