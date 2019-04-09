#' Split Newick Number String
#'
#' Takes a Newick number string and splits it into its' component parts.
#'
#' Prentheeses, commas, labels and semi-colon.
#'
#' @param NewickNumberString A vector of Newick Number Strings.
#'
#' @return A vector containing each component of a Newick number string.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get information for trees with N tips:
#' NewickNumberStringSplitter("(1,(2,(3)));")
#'
#' @export NewickNumberStringSplitter
NewickNumberStringSplitter <- function(NewickNumberString) {
  
  # Do basic string split first (may be incomplete as two or more digit numbers will be split in part):
  SplitString <- strsplit(NewickNumberString, split = "")[[1]]
  
  # Get position of numbers:
  NumberPositions <- grep("[:0-9:]", SplitString)
  
  # While any numbers run together (are split in part):
  while(any(diff(NumberPositions) == 1)) {
    
    # Find first split number:
    FirstSplitNumber <- NumberPositions[which(diff(NumberPositions) == 1)[1]]
    
    # Join number back togetehr and store:
    SplitString[FirstSplitNumber] <- paste(SplitString[FirstSplitNumber:(FirstSplitNumber + 1)], collapse = "")
    
    # Remove split part from string:
    SplitString <- SplitString[-(FirstSplitNumber + 1)]
    
    # Update number positions:
    NumberPositions <- grep("[:0-9:]", SplitString)
    
  }
  
  # Return split string with two or more digit numbers intact:
  return(SplitString)
  
}
