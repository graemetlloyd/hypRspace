#' Ladderize Newick Number String
#'
#' Takes as input Newick number string and returns the left ladderized version.
#'
#' @param NewickNumberString A split Newick number string.
#'
#' @return A vector of labelled Newick strings corresponding to unique trees.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get information for trees with N tips:
#' LadderizeNewickNumberString(NewickNumberStringSplitter("((4),(1,(2)));"))
#'
#' @export LadderizeNewickNumberString
LadderizeNewickNumberString <- function(NewickNumberString) {
  
  # MAKE FASTER/MORE ELEGANT!
  # REARRANGE SO DEALS WITH STRING NOT SPLIT STRING (OUTPUT TOO)
  
  # Find positions in string corresponding to numbers:
  NumberPositions <- grep("[:0-9:]", NewickNumberString)
  
  # Get values of those numbers:
  NumberValues <- as.numeric(NewickNumberString[NumberPositions])
  
  # Get total number of tips:
  NTips <- sum(NumberValues)
  
  # Make fake taxa labels:
  FakeTaxa <- MakeFakeLabels(NTips)
  
  # Create empty string matches vector for easy find and replace later:
  StringMatches <- vector(mode = "character")
  
  # For each number position:
  for(i in 1:length(NumberPositions)) {
    
    # Repalce number ith fake taxa:
    NewickNumberString[NumberPositions[i]] <- paste(FakeTaxa[1:NumberValues[i]], collapse = ",")
    
    # Add strings to vector (for easy find and replace later):
    StringMatches <- c(StringMatches, NewickNumberString[NumberPositions[i]])
    
    # Prne used taxa from fake taxa:
    FakeTaxa <- FakeTaxa[-(1:NumberValues[i])]
    
  }
  
  # Make into ape tree:
  tree <- ape::read.tree(text = paste(NewickNumberString, collapse = ""))
  
  # Get ladderized Newick:
  tree <- ape::write.tree(ape::ladderize(tree, right = FALSE))
  
  # Replace tips with numbers:
  for(i in StringMatches) tree <- gsub(i, as.character(length(strsplit(i, split = ",")[[1]])), tree)
  
  # Return in same split string format as input:
  return(strsplit(tree, split = "")[[1]])
  
}
