#' Get Unique Newick Number Strings
#'
#' Takes as input a vector of Newick number strings and returns only those pertaining to unique treeshapes.
#'
#' Principle of free rotation.
#'
#' @param NewickNumberStrings A vector of labelled Newick strings.
#'
#' @return A vector of Newick number strings corresponding to unique treeshapes.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get information for trees with N tips:
#' UniqueNewickNumberStrings(c("(1,(2),(2));", "(1,(2),(2));"))
#'
#' @export UniqueNewickNumberStrings
UniqueNewickNumberStrings <- function(NewickNumberStrings) {
  
  # Subfunction to add same label ("A") to each tip:
  AddSameLabelToAllTipsInNewickNumberString <- function(NewickNumberString) {
    
    # Turn Newick number into vector of characters:
    SplitString <- NewickNumberStringSplitter(NewickNumberString)
    
    # Get number positions:
    NumberPositions <- grep("[:0-9:]", SplitString)
    
    # Get number values:
    NumberValues <- as.numeric(SplitString[NumberPositions])
    
    # Replace numbers with N tips labelled A and separated by commas:
    SplitString[NumberPositions] <- unlist(lapply(as.list(NumberValues), function(x) paste(rep("A", times = x), collapse = ",")))
    
    # Recombine split string into single (now labelled) Newick string:
    return(paste(SplitString, collapse = ""))
    
  }
  
  # Add same label to all tips in each Newick number string:
  NewickStrings <- unlist(lapply(as.list(NewickNumberStrings), AddSameLabelToAllTipsInNewickNumberString))
  
  # Convert trees to an ape formatted multiPhylo object:
  Trees <- ape::read.tree(text = NewickStrings)
  
  # If a phylo object (single tree):
  if(class(Trees) == "phylo") {
    
    # Just return unmodiifed Newick string:
    return(NewickNumberStrings)
    
  # If a multiPhylo object:
  } else {
    
    # Add unit branch lengths to trees:
    Trees <- lapply(Trees, function(x) {x$edge.length <- rep(1, times = nrow(x$edge)); x})
    
    # Subfunction to convert an individual tree to a sorted string of it's splits:
    ConvertTreeToSortedSplitString <- function(Tree) {
      
      # Get number of tree tips:
      NTips <- ape::Ntip(Tree)
      
      # Get number of internal nodes:
      NInternalNodes <- ape::Nnode(Tree)
      
      # Get Internal node numbers:
      InternalNodeNumbers <- (NTips + 1):(NTips + NInternalNodes)
      
      # Add root age:
      Tree$root.time <- max(diag(ape::vcv(Tree)))
      
      # Get node ages:
      NodeAges <- unname(Claddis::GetNodeAges(Tree))
      
      # Return sorted string of tip ages and splits with appended node ages:
      return(paste(paste(sort(NodeAges[1:ape::Ntip(Tree)]), collapse = "%%"), paste(sort(paste(NodeAges[InternalNodeNumbers], unlist(lapply(as.list(InternalNodeNumbers), function(x) paste(Tree$tip.label[strap::FindDescendants(n = x, tree = Tree)], collapse = ""))), sep = "")), collapse = "%%"), sep = ""))
      
    }
    
    # Get sorted sorted string of tip ages and splits with appended node ages:
    TreesAsSortedSplitStrings <- unlist(lapply(Trees, ConvertTreeToSortedSplitString))
    
    # Add to each string the numbers from each Newick number string (creates truly unique unlabelled tree strings):
    TreesAsSortedSplitStrings <- paste(TreesAsSortedSplitStrings, unlist(lapply(lapply(as.list(NewickNumberStrings), NewickNumberStringSplitter), function(x) paste(sort(as.numeric(x[grep("[:0-9:]", x)])), collapse = "%%"))), sep = "")
    
    # Return tree strings with duplicates removed:
    return(NewickNumberStrings[!duplicated(TreesAsSortedSplitStrings)])
    
  }
  
}
