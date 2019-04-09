#' All Partition Combinations
#'
#' Given a set of partition sizes and labels produces all possible combinations.
#'
#' @param PartitionSizes A numeric vector of partition sizes.
#' @param Labels A character vector of labels.
#'
#' @return A matrix of all possible label combinations.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Get information for trees with N tips:
#' GetAllPartitionCombinations(c(1,2,3), LETTERS[1:6])
#'
#' @export GetAllPartitionCombinations
GetAllPartitionCombinations <- function(PartitionSizes, Labels) {
  
  # Create starting labels matrix:
  LabelsMatrix <- t(utils::combn(x = Labels, m = PartitionSizes[1]))
  
  # As long as there is more than one clade (i.e., it isn't the star tree):
  if(length(PartitionSizes) > 1) {
    
    # For each subsequent partition size in sequence:
    for(i in PartitionSizes[2:length(PartitionSizes)]) {
      
      # Increase abels matrix:
      LabelsMatrix <- do.call(rbind, lapply(apply(LabelsMatrix, 1, list), function(y) {CurrentLabels <- setdiff(Labels, unlist(y)); NewCombos <- t(combn(x = CurrentLabels, m = i)); cbind(do.call(rbind, rep(y, times = nrow(NewCombos))), NewCombos)}))
      
    }
    
  }
  
  # Return labels matrix:
  return(LabelsMatrix)
  
}
