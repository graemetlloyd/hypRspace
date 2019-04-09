#' Generate All Labelled Newicks from Newick Number Strings
#'
#' Generates all possible labelled newick strings from a set of Newick number strings.
#'
#' @param NewickNumberStrings A vector of Newick number strings.
#' @param TipLabels A character vector of tip labels to use.
#'
#' @return A list, one for each input Newick number string, each containing all possible labelled Newicks.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get information for trees with N tips:
#' AllLabelledNewicksFromNewickNumberStrings(c("(1,(2,(3)));", "(2,(1,(3)));"))
#'
#' @export AllLabelledNewicksFromNewickNumberStrings
AllLabelledNewicksFromNewickNumberStrings <- function(NewickNumberStrings, TipLabels = NULL) {
  
  # Get tip counts for each string:
  TipCounts <- unlist(lapply(as.list(NewickNumberStrings), function(x) {y <- NewickNumberStringSplitter(x); sum(as.numeric(y[grep("[:0-9:]", y)]))}))
  
  # If tip counts vary stop and warn user:
  if(length(unique(TipCounts)) > 1) stop("Tip counts are not consistent across Newick number strings.")
  
  # Set N as tip count:
  N <- unique(TipCounts)
  
  # Convert Newick strings to vectors of numbers:
  StringsAsNumericVectors <- lapply(strsplit(gsub("\\(|)|;", "", NewickNumberStrings), split = ","), as.numeric)
  
  # If no labels are provided then make fake ones:
  if(is.null(TipLabels)) TipLabels <- MakeFakeLabels(N)
  
  # Subfunction to build Newick strings from a Newick Number string and associated matrix of label combinations):
  BuildNewickStrings <- function(NewickNumberString, LabelsMatrix) {
    
    # Set up labelled newick strings by repeating Newick number string required number of times:
    LabelledNewickStrings <- rep(NewickNumberString, times = nrow(LabelsMatrix))
    
    # Create split strings in list format ready for mapply:
    SplitStrings <- lapply(as.list(LabelledNewickStrings), NewickStringSplitter)
    
    # Create split string labels as list format ready for mapply:
    StringLabels <- lapply(apply(LabelsMatrix, 1, list), unlist)
    
    # Subfunction to add labels to Newick number string:
    LabelNewickNumberString <- function(SplitString, StringLabels) {
      
      # Find positions of numbers in split string:
      NumberPositions <- grep("[:0-9:]", SplitString)
      
      # Find number values for split string (N labels required):
      NumberValues <- as.numeric(SplitString[NumberPositions])
      
      # For each number positions (CAN PROBABLY REPLACE THIS WITH A SINGLE LINE SOMEHOW):
      for(i in 1:length(NumberPositions)) {
        
        # Paste string labels over ith number position:
        SplitString[NumberPositions[i]] <- paste(StringLabels[1:as.numeric(NumberValues[i])], collapse = ",")
        
        # Remove used string labels from pool:
        StringLabels <- StringLabels[-c(1:as.numeric(NumberValues[i]))]
        
      }
      
      # Convert split string back into a single string:
      paste(SplitString, collapse = "")
      
    }
    
    # Generate all labelled Newick strings from Newick number string and labels:
    LabelledNewickStrings <- mapply(LabelNewickNumberString, SplitString = SplitStrings, StringLabels = StringLabels)
    
    # Prune duplicated Newick strings:
    LabelledNewickStrings <- UniqueNewickStrings(LabelledNewickStrings)
    
    # Return labelled Newick strings:
    return(LabelledNewickStrings)
    
  }
  
  # Generate and return all labelled trees:
  mapply(BuildNewickStrings, NewickNumberString = as.list(NewickNumberStrings), LabelsMatrix = lapply(StringsAsNumericVectors, GetAllPartitionCombinations, Labels = TipLabels), SIMPLIFY = FALSE)
  
}
