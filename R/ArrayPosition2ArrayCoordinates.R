#' Get Array Coordinates From Array Position
#'
#' Simply returns the coordinates for the Nth value in an array.
#'
#' Depending on different criteria N tips will produce a specific number of trees.
#'
#' @param Position The number of the Nth value in the array (an integer).
#' @param ArrayDimensions A vector of integers giving the dimensions of the array.
#'
#' @return A vector the same length as \code{ArrayDimensions} giving the coordinates.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @keywords array, coordinates
#'
#' @examples
#'
#' # Get the coordinates of the 17th value in a
#' # 2 by 4 by 6 array:
#' ArrayPosition2ArrayCoordinates(17, c(2,4,6))
#'
#' @export ArrayPosition2ArrayCoordinates
ArrayPosition2ArrayCoordinates <- function(Position, ArrayDimensions) as.vector(which(array(1:prod(ArrayDimensions), dim = ArrayDimensions) == Position, arr.ind = TRUE))
