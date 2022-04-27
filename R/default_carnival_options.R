#' Setting Default CARNIVAL Options
#'
#' Returns the default CARNIVAL options as a list.  You can modify the elements
#' of the list and then use it as an argument in \code{\link{run_phonemes}}.
#' If you choose CPLEX or CBC, you must modify then the solverPath field and point to
#' the CPLEX/CBC executable (See Details).
#'
#' @details PHONEMeS is dependent on CARNIVAL for exhibiting the signalling pathway optimisation
#' which requires the interactive version of IBM Cplex, Gurobi or CBC-COIN solver
#'  as the network optimiser. The IBM ILOG Cplex is freely available through
#' Academic Initiative \href{https://www.ibm.com/products/ilog-cplex-optimization-studio}{here}.
#' Gurobi license is also free for academics, request a license following instructions
#' \href{https://www.gurobi.com/downloads/end-user-license-agreement-academic/}{here}.
#' The  \href{https://projects.coin-or.org/Cbc}{CBC} solver is open source and freely
#' available for any user, but has a significantly lower performance than CPLEX or
#' Gurobi. Obtain CBC executable directly usable for cosmos
#'\href{https://ampl.com/products/solvers/open-source/#cbc}{here}. Alternatively for
#' small networks, users can rely on the freely available
#'  \href{https://cran.r-project.org/web/packages/lpSolve/index.html}{lpSolve R-package},
#'  which is automatically installed with the package.
#'
#' @param solver one of `cplex` (recommended, but require 3rd party tool), `cbc` (also require 3rd party tool) or `lpSolve` (only for small networks)
#' @return returns a list with all possible options implemented in CARNIVAL
#' @examples
#' # load and change default options:
#' my_options = default_carnival_options(solver = "cplex")
#'
#' my_options$solverPath = "/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex"
#' my_options$threads = 2
#' my_options$timelimit = 3600*15
#' @export
#'
default_carnival_options <- function(solver = NULL){

  if(is.null(solver)){
    stop("please call default_carnival_options(solver = 'cplex') with a
        specific solver argument. Valid solvers are: 'lpSolve','cplex' or 'cbc'.")
  }
  solver_options <- c("lpSolve","cplex","cbc")
  solver <- match.arg(solver,choices = solver_options)

  if(solver == "lpSolve"){
    opts <- CARNIVAL::defaultLpSolveCarnivalOptions()
  }else if(solver == "cplex"){
    opts <- CARNIVAL::defaultCplexCarnivalOptions()
  }else if(solver == "cbc"){
    opts <- CARNIVAL::defaultCbcSolveCarnivalOptions()
  }
  opts$keepLPFiles <- FALSE

  # for older Carnival versions
  opts$outputFolder <- getwd()
  opts$workdir <- getwd()

  return(opts)
}



#' Check CARNIVAL Options
#'
#' Checks the list of input options for CARNIVAL for completeness.
#'
#' @param opts An object of type \dQuote{\code{list}} defining the run parameters
#'   for CARNIVAL  Use the \code{\link{default_carnival_options}} function to
#'   create a list with default parameter settings. If cplex or cbc are chosen
#'   as the solver, the parameter solverPath needs to be supplied (not
#'   automatically added by \code{default_carnival_options()}).
#'
#' @seealso \code{\link{default_carnival_options}},
#'   \code{\link[PHONEMeS]{run_phonemes}}
#' @noRd
check_carnival_options <- function(opts){


  if(!is.list(opts)) stop("CARNIVAL options should be a list")
  req_names <- c("solver")


  if(!all(req_names %in% names(opts))){
    stop("CARNIVAL options should contain all options.
             Start by calling PHONEMeS::default_carnival_options() with your solver specified. ")
  }


  if(!(opts$solver %in% c("cplex","cbc","lpSolve"))) stop("current version supports only the CPLEX or cbc solver. lpSolve is also available for test runs.")
  if(opts$solver == "cbc") print("We recommend the users to use CPLEX if possible.")
  if(opts$solver == "lpSolve") print("lpSolve does not scale well with large PKNs. This solver is mainly for testing purposes. To run CARNIVAL for PHONEMeS, we recommend using cplex, or cbc solvers.")
  if(is.null(opts$solverPath) & opts$solver != "lpSolve") stop("path to CPLEX or cbc solver must be provided")

  return('List of CARNIVAL options is complete.')


}
