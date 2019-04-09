#' Split Newick String
#'
#' Takes a Newick string and splits it into its' component parts.
#'
#' Parentheeses, commas, labels and semi-colon.
#'
#' @param NewickString A vector of Newick strings.
#'
#' @return A vector containing each component of a Newick string.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get information for trees with N tips:
#' NewickStringSplitter("(Aa,(Bb,Cc,(Dd,Ee,Ff)));")
#'
#' @export NewickStringSplitter
NewickStringSplitter <- function(NewickString) {
  
  # Do basic string split first (may be incomplete as two or more digit numbers will be split in part):
  SplitString <- strsplit(NewickString, split = "")[[1]]
  
  # Get positions of labels:
  LabelPositions <- which(!apply(cbind(SplitString == "(", SplitString == ")", SplitString == ",", SplitString == ";"), 1, any))
  
  # While any labels run together (are split in part):
  while(any(diff(LabelPositions) == 1)) {
    
    # Find first split label:
    FirstSplitNumber <- LabelPositions[which(diff(LabelPositions) == 1)[1]]
    
    # Join label back togetehr and store:
    SplitString[FirstSplitNumber] <- paste(SplitString[FirstSplitNumber:(FirstSplitNumber + 1)], collapse = "")
    
    # Remove split part from string:
    SplitString <- SplitString[-(FirstSplitNumber + 1)]
    
    # Update label positions:
    LabelPositions <- which(!apply(cbind(SplitString == "(", SplitString == ")", SplitString == ",", SplitString == ";"), 1, any))
    
  }
  
  # Return split string with labels intact:
  return(SplitString)
  
}
