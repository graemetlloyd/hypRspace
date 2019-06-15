#' Generate All Unlabelled Trees Through Partitions
#'
#' Generates all unlabelled rooted multifurcating trees of N tips.
#'
#' Function works using notion in Felsenstein (2004) of applying partitions (combinatorics) to generate all possible sets of splits.
#'
#' Note that currently this function fails to generate all possible trees past N = 7 so it is not recommended for use.
#'
#' @param N The required number of tips.
#'
#' @return A list of Newick number strings partitioned by number of internal nodes (1 to N - 1).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @keywords combinatorics, phylogeny, trees
#'
#' @references
#'
#' Felsenstein 2004
#'
#' @examples
#'
#' # Generate all trees of five tips:
#' AllNewickNumbersByPartition(N = 5)
#'
#' @export AllNewickNumbersByPartition
AllNewickNumbersByPartition <- function(N) {
  
  # Subfunction to generate phylogenetic partitions:
  BuildPhyloPartitions <- function(N) {
    
    # Get raw partitions of N:
    RawPartitions <- lapply(lapply(apply(partitions::parts(N), 2, list), function(x) {x <- unlist(x); x[x != 0]}), rev)
    
    # Find any partitions that are all ones (effectively same as star phylogeny and hence needs to be removed):
    AllOnes <- which(unlist(lapply(RawPartitions, function(x) all(x == 1))))
    
    # Remove all one partitions:
    if(length(AllOnes) > 0) RawPartitions <- RawPartitions[-AllOnes]
    
    # Build initial partial Newick number strings:
    PhyloPartitions <- lapply(lapply(RawPartitions, function(x) c(if(any(x == 1)) as.character(sum(x == 1)), paste("(", x[x > 1], ")", sep = ""))), paste, collapse = ",")
    
    # Output phylo partitions:
    return(PhyloPartitions)
    
  }
  
  # Only proceed if at least two leaves (otherwise it is only star trees - covered below):
  if(N > 2) {
    
    # Generate all aprtial Newick number strings (excludes "outer" root):
    PartialNewickNumberStrings <- BuildPhyloPartitions(N)
    
    # Add first (i.e., star) phylogeny to final Newick number strings:
    FinalNewickNumberStrings <- paste(PartialNewickNumberStrings[[1]], ";", sep = "")
    
    # Remove star phylogeny and store partial newick number strings as a vector:
    PartialNewickNumberStrings <- unlist(PartialNewickNumberStrings[2:length(PartialNewickNumberStrings)])
    
    # While there are still partial Newick strings:
    while(length(PartialNewickNumberStrings) > 0) {
      
      # Get split version of first partial Newick number string:
      SplitString <- NewickNumberStringSplitter(PartialNewickNumberStrings[1])
      
      # Get values of each number in the partial Newick number string (i.e., to check whether more partitions are possible):
      NumberValues <- as.numeric(SplitString[grep("[:0-9:]", SplitString)])
      
      # If any partitions are greater than 2 in size (more partitions possible):
      if(any(NumberValues > 2)) {
        
        # Create paritions of each unique number required.
        AllPhyloPartitionsForEachNumberValueOverTwo <- lapply(as.list(unique(NumberValues[NumberValues > 2])), function(x) unlist(BuildPhyloPartitions(x)))
        
        # If there is only one position where partitions need to be permuted:
        if(sum(NumberValues > 2) == 1) {
          
          # Find position of number that needs partitioning:
          NumberPosition <- which(SplitString == as.character(NumberValues[NumberValues > 2]))
          
          # Build new strings by inserting each partition into the number position:
          NewStrings <- c(paste("(", PartialNewickNumberStrings[1], ");", sep = ""), paste("(", unlist(lapply(mapply(function(x, y) {x[NumberPosition] <- y; x}, x = rep(list(SplitString), length(AllPhyloPartitionsForEachNumberValueOverTwo[[1]]) - 1), y = AllPhyloPartitionsForEachNumberValueOverTwo[[1]][-1], SIMPLIFY = FALSE), paste, collapse = "")), ");", sep = ""))
          
        # If two or more spots where new partitions need to be permuted:
        } else {
          
          # Add names to AllPhyloPartitionsForEachNumberValueOverTwo:
          names(AllPhyloPartitionsForEachNumberValueOverTwo) <- unique(NumberValues[NumberValues > 2])
          
          # Build out all new partitions:
          NewPartitions <- lapply(as.list(NumberValues[NumberValues > 2]), function(x) AllPhyloPartitionsForEachNumberValueOverTwo[[which(names(AllPhyloPartitionsForEachNumberValueOverTwo) == x)]])
          
          # Build array of string with each "axis" representing a number that needs partitions exploring:
          StringArray <- array(rep(list(SplitString), prod(unlist(lapply(NewPartitions, length)))), dim = unlist(lapply(NewPartitions, length)), dimnames = NewPartitions)
          
          # Find position of numbers to partition in split string:
          NumberPositions <- grep("[:0-9:]", SplitString)[NumberValues > 2]
          
          # Strip out parentheses in first name (star phylogeny) to avoid these duplicating later (i.e., want "(5)" not "((5))"):
          for(i in 1:length(dim(StringArray))) dimnames(StringArray)[[i]][1] <- gsub("\\(|\\)", "", dimnames(StringArray)[[i]][1])
          
          # Get array dimensions:
          NArrayDimensions <- length(dim(StringArray))
          
          # For each element in the string array:
          for(i in 1:length(StringArray)) {
            
            # Find array coordinates of ith element:
            ArrayCoordinates <- ArrayPosition2ArrayCoordinates(Position = i, ArrayDimensions = dim(StringArray))
            
            # Craete empty partitions to substitute vector:
            PartitionsToSubstitute <- vector(mode = "character")
            
            # Populate partitions to substitute vector:
            for(j in 1:NArrayDimensions) PartitionsToSubstitute <- c(PartitionsToSubstitute, dimnames(StringArray)[[j]][ArrayCoordinates[j]])
            
            # Insert new partitions in correct part of split string:
            StringArray[i][[1]][NumberPositions] <- PartitionsToSubstitute
            
            # Format as complete string and add abck to array:
            StringArray[i] <- paste(StringArray[i][[1]], collapse = "")
            
          }
          
          # Build new strings list by adding root node and end line:
          NewStrings <- paste("(", unlist(StringArray), ");", sep = "")
          
          # Remove any accidentally generated duplicate strings:
          NewStrings <- UniqueNewickNumberStrings(NewStrings)
          
        }
        
        # Add fully formatted Newick number string to final list:
        FinalNewickNumberStrings <- c(FinalNewickNumberStrings, setdiff(NewStrings, FinalNewickNumberStrings))
        
        # Find any strings that need further partitioning:
        StringsToKeepPartitioning <- which(unlist(lapply(as.list(NewStrings), function(x) max(as.numeric(NewickNumberStringSplitter(x)[grep("[:0-9:]", NewickNumberStringSplitter(x))])))) > 2)
        
        # If strings to keep partitioning exist add them to list:
        if(length(StringsToKeepPartitioning) > 0) PartialNewickNumberStrings <- c(PartialNewickNumberStrings, setdiff(unlist(lapply(as.list(NewStrings[StringsToKeepPartitioning]), function(x) paste(NewickNumberStringSplitter(x)[2:(length(NewickNumberStringSplitter(x)) - 2)], collapse = ""))), PartialNewickNumberStrings))
        
        # Remove partial Newick number string from the list (allows move to next string):
        PartialNewickNumberStrings <- PartialNewickNumberStrings[-1]
        
      # If string cannot be partitioned further (all remaining partitions are 2 or less in size):
      } else {
        
        # Add fully formatted Newick number string to final list:
        FinalNewickNumberStrings <- c(FinalNewickNumberStrings, paste("(", PartialNewickNumberStrings[1], ");", sep = ""))
        
        # Remove partial Newick number string from the list (allows move to next string):
        PartialNewickNumberStrings <- PartialNewickNumberStrings[-1]
        
      }
      
    }
    
    # Collapse to just unique strings:
    FinalNewickNumberStrings <- UniqueNewickNumberStrings(FinalNewickNumberStrings)
    
    # Get number of internal nodes in each tree:
    NInternalNodes <- unlist(lapply(as.list(FinalNewickNumberStrings), function(x) sum(NewickNumberStringSplitter(x) == "(")))
    
    # Reorganise output into list of trees with each node count from 1 to N - 1:
    FinalNewickNumberStrings <- lapply(as.list(unique(NInternalNodes)), function(x) FinalNewickNumberStrings[x == NInternalNodes])
    
  # If N is two or less:
  } else {
    
    # Build one leaf Newick number string:
    if(N == 1) FinalNewickNumberStrings <- list("(1);")
    
    # Build two leaf Newick number string:
    if(N == 2) FinalNewickNumberStrings <- list("(2);")
    
  }
  
  # Output final Newick number strings:
  return(FinalNewickNumberStrings)
  
}

