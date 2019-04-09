#' Get Unique Newick Strings
#'
#' Takes as input a vector of labelled newick strings and reqturns only those pertaining to unique trees.
#'
#' Principle of free rotation.
#'
#' @param NewickStrings A vector of labelled Newick strings.
#'
#' @return A vector of labelled Newick strings corresponding to unique trees.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get information for trees with N tips:
#' UniqueNewickStrings(c("(A,(B,C),(D,E));", "(A,(E,D),(C,B));"))
#'
#' @export UniqueNewickStrings
UniqueNewickStrings <- function(NewickStrings) {
  
  # Convert trees to an ape formatted multiPhylo object:
  Trees <- ape::read.tree(text = NewickStrings)
  
  # If a phylo object (single tree):
  if(class(Trees) == "phylo") {
    
    # Just return unmodiifed Newick string:
    return(NewickStrings)
    
  # If a multiPhylo object:
  } else {
    
    # Subfunction to convert an individual tree to a sorted string of it's splits:
    ConvertTreeToSortedSplitString <- function(Tree) {
      
      # Get number of tree tips:
      NTips <- ape::Ntip(Tree)
      
      # Get number of internal nodes:
      NInternalNodes <- ape::Nnode(Tree)
      
      # Get Internal node numbers:
      InternalNodeNumbers <- (NTips + 1):(NTips + NInternalNodes)
      
      # Return soretd string of splits:
      return(paste(sort(unlist(lapply(as.list(InternalNodeNumbers), function(x) paste(sort(Tree$tip.label[strap::FindDescendants(n = x, tree = Tree)]), collapse = ",")))), collapse = "%%"))
      
    }
    
    # Get sorted tree string of splits for each tree (duplicates here are duplicated trees in free rotation):
    TreesAsSortedSplitStrings <- unlist(lapply(Trees, ConvertTreeToSortedSplitString))
    
    # Return tree strings with duplicates removed:
    return(NewickStrings[!duplicated(TreesAsSortedSplitStrings)])
    
  }
  
}
