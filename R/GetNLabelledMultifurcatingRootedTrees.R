#' Get N Labelled Multifucrating Rooted Trees
#'
#' Gets the total number of labelled multifurcating rooted trees for N tips.
#'
#' Uses equation 3.1 (and illustrated in Figure 3.6) of Felsenstein (2004).
#'
#' Note that direct counts are possible only up to N = 145, beyond that R will return only "Inf".
#'
#' @param N The number of tips.
#'
#' @return A list containg two items, TotalTrees (the total number of trees) and TreesByNodeCount (a vector of trees by node count from 1 to N - 1).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @keywords combinatorics, phylogeny, trees
#'
#' @see also
#'
#' The function ape::howmanytrees attempts something similar, but has more options. Whereas phangorn::allTrees actually generates trees, but only the bifurcating ones and only up to 10 tips.
#'
#' @references
#'
#' Felsenstein 2004
#'
#' @examples
#'
#' # Get numbers of trees:
#' GetNLabelledMultifurcatingRootedTrees(8)
#'
#' # Get just the total number of trees for 1:10 tips
#' unlist(lapply(as.list(1:10), function(x)
#'   GetNLabelledMultifurcatingRootedTrees(x)$TotalTrees))
#'
#' @export GetNLabelledMultifurcatingRootedTrees
GetNLabelledMultifurcatingRootedTrees <- function(N) {
  
  # Set default (starting) values for N = 1, or N = 2):
  TreesByNodeCount <- TotalTrees <- 1
  
  # Only eed to modify if requested N is greater than 2):
  if(N > 2) {
    
    # Iteratively add trees using equation 3.1 of Felsenstein (2004):
    for(i in 3:N) TreesByNodeCount <- c(TreesByNodeCount * (1:(i - 2)), 0) + c(0, TreesByNodeCount * (i:(i + length(TreesByNodeCount) - 1)))
    
    # Update total trees:
    TotalTrees <- sum(TreesByNodeCount)
    
  }
  
  # Compile output in a list:
  output <- list(TotalTrees, TreesByNodeCount)
  
  # Add anems to output:
  names(output) <- c("TotalTrees", "TreesByNodeCount")
  
  # Return output:
  return(output)
  
}

# Make plot of all total tree solutions that can be directly calculated:
#plot(1:145, unlist(lapply(as.list(1:145), function(x) GetNLabelledMultifurcatingRootedTrees(x)$TotalTrees)), xlab = "Leaf count", ylab = "Tree count", log = "y", type = "l")
