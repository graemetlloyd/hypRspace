#' Make fake labels
#'
#' Generates a set of fake labels (character strings).
#'
#' @param N The number of labels required.
#'
#' @return A vector of labels (character strings).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get information for trees with N tips:
#' MakeFakeLabels(N = 100)
#'
#' @export MakeFakeLabels
MakeFakeLabels <- function(N) {
  
  # Subfunction to generate fake labels:
  FakeLabelRecursor <- function(y) unlist(lapply(lapply(as.list(y), rep, length.out = 26), function(y) paste(y, LETTERS, sep = "")))
  
  # If more tips than letters:
  if(N > 26) {
    
    # Generate fake labels:
    FakeLabels <- FakeLabelRecursor(LETTERS)
    
    # Keep generating labels until there are at least N:
    while(N > length(FakeLabels)) FakeLabels <- FakeLabelRecursor(FakeLabels)
    
  # If no more taxa required than there are letters:
  } else {
    
    # Set fake taxa as letters:
    FakeLabels <- LETTERS
    
  }
  
  # Reduce fake labels to just N:
  FakeLabels <- FakeLabels[1:N]
  
  # Return fake taxa:
  return(FakeLabels)
  
}
