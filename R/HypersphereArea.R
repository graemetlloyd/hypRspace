#' Get hypersphere area
#'
#' Calculates the dimensional equivalent of area for a hypersphere.
#'
#' Given a specified radius and number of dimensions calculates the dimensional equivalent of area for a hypersphere.
#'
#' @param Radius The radius of the sphere.
#' @param NDimensions The number of dimensions of the sphere.
#'
#' @return The surface area (or dimensional equivalent of area).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @keywords hypersphere
#'
#' @examples
#'
#' # Get surface area of the regular unit sphere:
#' HypersphereArea <- function(Radius = 1, NDimension = 3)
#'
#' @export HypersphereArea
HypersphereArea <- function(Radius, NDimensions) {
  
  # Subfunction to calculate the double factorial (peel off into own funtion maybe?):
  DoubleFactorial <- function(N) {
    
    # If N is -1 return 1:
    if(N == -1) x <- 1
    
    # If N is 0 return 1:
    if(N == 0) x <- 1
    
    # If N is odd return double factorial:
    if(N %% 2 == 1) x <- prod(which(1:N %% 2 == 1))
    
    # If N is even return double factorial:
    if(N %% 2 == 0) x <- prod(which(1:N %% 2 == 0))
    
    # Return double factorial:
    return(x)
    
  }
  
  # If odd number of dimensions:
  if(NDimensions %% 2 == 1) {
    
    # Calculate area of unit sphere:
    S_n <- ((2 ^ ((NDimensions + 1) / 2)) * (pi ^ ((NDimensions - 1) / 2))) / DoubleFactorial(NDimensions - 2)
    
  }
  
  # If even number of dimensions:
  if(NDimensions %% 2 == 0) {
    
    # Calculate area of unit sphere:
    S_n <- (2 * (pi ^ (NDimensions / 2))) / factorial((0.5 * NDimensions) - 1)

  }
  
  # Return area scaled by radius to power of N dimensions minus one:
  return(S_n * (Radius ^ (NDimensions - 1)))
  
}
