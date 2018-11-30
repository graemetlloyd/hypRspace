#' Get N Trees
#'
#' Gets the total number of trees for N tips.
#'
#' Depeding on different criteria N tips will produce a specific number of trees.
#'
#' @param N The number of tips.
#'
#' @return {The number of trees.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @keywords trees
#'
#' @examples
#'
#' # Create some example coordinates:
#' N <- 3
#'
#' @export GetNTrees
GetNTrees <- function(N) {
  
  # N bifurcating
  # N any furcation
  # N possible partitions
  # N unique tree shapes (unlabelled trees)
  
  # Function to calculate total number of (labelled) bipartitions for N taxa (excludes root and tips, i.e., m = 1 or m = N):
  NBipartitions <- function(N) return(sum(unlist(lapply(lapply(as.list((N - 1):2), function(x) combn(x = N, m = x)), ncol))))
  
  
  
  #treeshapes is just nested combn calls...
  #factorial((2 * factorial(N)) - 3) / ((2 ^ (N - 2)) * factorial(N - 2))
  
  return(NBipartitions)
  
  # EXACT VERSUS ESTIMATED IN OUTPUT
  # DO NOT REALLY NEED TO CLALCULATE AS CAN JUST HAVE A TABLE FOR ALL DIRECTLY ESTIMABLE VALUES

# Estimate as best as possible for larger numbers (e.g., to nearest 10 ^ N). (Also, work out at which level regular ones fail.)

}
