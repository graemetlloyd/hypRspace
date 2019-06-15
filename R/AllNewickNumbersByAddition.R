#' Generate All Unlabelled Trees Through Addition
#'
#' Generates all unlabelled rooted multifurcating trees of N tips.
#'
#' Function works using notion in Felsenstein (2004) or simply starting with the smallest possible tree and additively generating trees up to N tips by adding a new tip in every possible place.
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
#' AllNewickNumbersByAddition(N = 5)
#'
#' @export AllNewickNumbersByAddition
AllNewickNumbersByAddition <- function(N) {
  
  ################################################################
  #         VARIOUS UNDOCUMENTED SUBFUNCTIONS BEGIN HERE         #
  ################################################################
  
  # Subfunction to find any "comma" numbers (unique restrictions apply):
  CommaNumberFinder <- function(x) {
    
    # Get positions before a comma in string:
    BeforeComma <- which(x == ",") - 1
    
    # Isolate just commas that follow a number (i.e., those that do not follow a close parentheses, the only other option):
    CommaNumbers <- BeforeComma[!x[BeforeComma] == ")"]
    
    # Output positions of numbers followed by a comma (if any):
    return(CommaNumbers)
    
  }
  
  # Subfunction to find any "parenthesis" numbers (unique restrictions apply):
  ParentheticNumberFinder <- function(x) {
    
    # Get positions before a close parenthesis in string:
    BeforeClosingParenthesis <- which(x == ")") - 1
    
    # Isolate just parentheses that follow a number (i.e., those that do not follow a close parentheses, the only other option):
    ParenthesisNumbers <- BeforeClosingParenthesis[!x[BeforeClosingParenthesis] == ")"]
    
    # Output positions of numbers inside parentheses:
    return(ParenthesisNumbers)
    
  }
  
  # Subfunction to get parenthesis matches:
  GetParenthesisMatches <- function(x) {
    
    # Create empty vectors to store nestedness, starting and ending parenthesis positions and temporary holding position:
    Nestedness <- StartPositions <- EndPositions <- HoldingPosition <- vector(mode = "numeric")
    
    # For each value from left to right in the Newick string:
    for(CurrentPosition in 1:(length(x) - 1)) {
      
      # If an opening parenthesis then add to the holding stack:
      if(x[CurrentPosition] == "(") HoldingPosition <- c(HoldingPosition, CurrentPosition)
      
      # If a closing parenthesis:
      if(x[CurrentPosition] == ")") {
        
        # Add to end positions:
        EndPositions <- c(EndPositions, CurrentPosition)
        
        # Pull off most recent opening parenthesis and add to start positions:
        StartPositions <- c(StartPositions, HoldingPosition[length(HoldingPosition)])
        
        # Record nestedness of parenthesis pair:
        Nestedness <- c(Nestedness, length(HoldingPosition) - 1)
        
        # Remove that starting parenthesis from the stack:
        HoldingPosition <- HoldingPosition[-length(HoldingPosition)]
        
      }
      
    }
    
    # Return start and end points:
    return(cbind(StartPositions, EndPositions, Nestedness))
    
  }
  
  # Subfunction to add a new tip to a numeric Newick string at opening parenthesis position(s):
  AddAtParentheses <- function(x) return(unique(unlist(lapply(apply(GetParenthesisMatches(x), 1, as.list), function(y) {y <- unlist(y); x[y[1]] <- "(1,("; x[y[2]] <- "))"; paste(x, collapse = "")}))))
  
  # Subfunction to permute all possible tree strings:
  TreeStringPermuter <- function(x, mode) {
    
    # Create empty vector to store output tree strings:
    OutputTreeStrings <- vector(mode = "character")
    
    # If mode is adding a tip at a node (no new internal nodes):
    if(mode == "AddAtNode") {
      
      # Get positions of numbers inside parentheses:
      ParenthesisPositions <- ParentheticNumberFinder(x)
      
      # Permute all possible Newick-like strings where no new internal node is added:
      OutputTreeStrings <- c(OutputTreeStrings, unique(unlist(lapply(as.list(grep("[:0-9:]", x)), function(y) {x[y] <- as.character(as.numeric(x[y]) + 1); paste(x, collapse = "")}))))
      
      # If all numbers correspond to clades (then possible to add single taxon in polytomy with clades):
      if(length(setdiff(grep("[:0-9:]", x), ParenthesisPositions)) == 0) {
        
        # Add extra taxon in polytomy with clades:
        x[ParenthesisPositions[1] - 1] <- paste("1,", x[ParenthesisPositions[1] - 1], sep = "")
        
        # Add Newick string to output:
        OutputTreeStrings <- c(OutputTreeStrings, paste(x, collapse = ""))
        
      }
      
    }
    
    # If mode is adding a tip (creating a new node):
    if(mode == "AddAtEdge") {
      
      # Get positions of numbers followed by a comma (if any):
      CommaPostitions <- CommaNumberFinder(x)
      
      # Add new tip additions to output tree strings:
      OutputTreeStrings <- c(OutputTreeStrings, AddAtParentheses(x))
      
      # If there are comma positions:
      if(length(CommaPostitions) > 0) {
        
        # Take all comma position numbers and add one whilst placing them in parentheses:
        OutputTreeStrings <- c(OutputTreeStrings, unique(unlist(lapply(as.list(CommaPostitions), function(y) {x[y] <- paste("(", as.numeric(x[y]) + 1, ")", sep = ""); paste(x, collapse = "")}))))
        
      }
      
    }
    
    # Prune ladderized duplicates:
    OutputTreeStrings <- unique(unlist(lapply(as.list(OutputTreeStrings), function(y) paste(LadderizeNewickNumberString(strsplit(y, split = "")[[1]]), collapse = ""))))
    
    # Return new output tree strings as a vector of numeric Newicks:
    return(OutputTreeStrings)
    
  }
  
  ################################################################
  #          VARIOUS UNDOCUMENTED SUBFUNCTIONS END HERE          #
  ################################################################
  
  # Build start tree (Newick but with numbers of taxa only):
  TreeStrings <- list("(1);")
  
  # For each N from 2 to desired N:
  for(i in 2:N) {
    
    # Store previous string so can still call whilst modifying TreeStrings:
    PreviousTreeStrings <- TreeStrings
    
    # Set maximum number of internal nodes for N tips:
    NInternalNodes <- max(c(1, i - 1))
    
    # For new internal node (if required) begin by populating new level with tree(s) from lower level:
    if(i > 2) TreeStrings[[NInternalNodes]] <- TreeStrings[[(NInternalNodes - 1)]]
    
    # For each level of resolution (from one internal node to N - 1 internal nodes):
    for(j in 1:NInternalNodes) {
      
      # Special case of dealing with a single (root) node (i.e., the star phylogeny):
      if(j == 1) {
        
        # Update single node start tree to just N in parentheses:
        TreeStrings[[j]] <- paste("(", i, ");", sep = "")
        
      # If dealing with at least two internal nodes:
      } else {
        
        # Build out split strings (each character an item in a vector) from previous N at same level:
        SplitString <- lapply(lapply(TreeStrings[[j]], NewickNumberStringSplitter), unlist)
        
        # Case if no internal nodes need to be added:
        if(NInternalNodes == j) {
          
          # Permute all Newick-like strings where a new internal node is added:
          TreeStrings[[j]] <- unique(unlist(lapply(SplitString, TreeStringPermuter, mode = "AddAtEdge")))
          
        # Case if new internal node needs to be added (changes how addition works):
        } else {
          
          # Permute all Newick-like strings where a tip is added at an existing internal node:
          TreeStrings[[j]] <- unique(unlist(lapply(SplitString, TreeStringPermuter, mode = "AddAtNode")))
          
        }
        
        # Build out split strings (each character an item in a vector) from previous N at lower level:
        SplitString <- lapply(lapply(PreviousTreeStrings[[(j - 1)]], NewickNumberStringSplitter), unlist)
        
        # Add a new node to these and add to TreeStrings any new unique unlabelled Newicks:
        TreeStrings[[j]] <- unique(c(TreeStrings[[j]], unique(unlist(lapply(SplitString, TreeStringPermuter, mode = "AddAtEdge")))))
        
      }
      
    }
    
  }
  
  # Prune out any duplicate Newick number strings (avoid principle of free rotation issue):
  TreeStrings <- lapply(TreeStrings, UniqueNewickNumberStrings)
  
  # Return generated unlabelled topologies as "numeric" Newick strings:
  return(TreeStrings)
  
}

