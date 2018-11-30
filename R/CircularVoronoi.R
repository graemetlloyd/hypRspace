#' Circular Voronoi Plot
#'
#' Makes a circular Voronoi plot.
#'
#' Intended for visualising two-dimensional projections of hyperdimensional point clouds that approximate a hypersphere. Is effectively a circular extension of the Voronoi function in the \link{deldir} package.
#'
#' @param x The x coordinates of the points to plot.
#' @param y The y coordinates of the points to plot.
#' @param heat Numerical values that will determine the colour of each tile plotted (i.e., as a heat map).
#' @param EdgeFactor The proportional extension of the edge of the plotting circle (i.e., the default of 1.1 adds 10% to the dimatere of the point cloud).
#' @param NColours The number of (possible) colours to use in plotting (default is 100).
#'
#' @return A circular Voronoi plot with tiles coloured according to heat values.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @keywords Voronoi
#'
#' @examples
#'
#' # Create some random data:
#' N <- 100
#' Radius <- runif(n = N, min = 0, max = 2)
#' Angle <- sample(1:360, size = N, replace = TRUE)
#' x <- Radius * cos(Angle * (pi / 180))
#' y <- Radius * sin(Angle * (pi / 180))
#' heat <- sort(runif(n = N))
#'
#' # Create treespace tile plot:
#' CircularVoronoi(x, y, heat, EdgeFactor = 1.05)
#'
#' @export CircularVoronoi
CircularVoronoi <- function(x, y, heat, EdgeFactor = 1.1, NColours = 100) {
  
  # ADD HEATMAP LEGEND
  # ADD COLOUR PALETTE OPTIONS
  
  # Check edge factor is a valid value and stop and warn if not:
  if(EdgeFactor < 1) stop("EdgeFactor must be equal to or greater than 1")
  
  # Set total number of values:
  N <- length(x)
  
  # Get Delaunay Triangulation from x and y:
  DelaunayForVoronoi <- deldir::deldir(x, y)$dirsgs
  
  # Add bearings (from origin) to each coordinate as last two columns:
  DelaunayForVoronoi <- cbind(DelaunayForVoronoi, ((180 / pi) * atan2(as.numeric(DelaunayForVoronoi[, "y1"]), as.numeric(DelaunayForVoronoi[, "x1"]))) %% 360, ((180 / pi) * atan2(as.numeric(DelaunayForVoronoi[, "y2"]), as.numeric(DelaunayForVoronoi[, "x2"]))) %% 360)
  
  # Update column names accordingly:
  colnames(DelaunayForVoronoi)[(ncol(DelaunayForVoronoi) - 1):ncol(DelaunayForVoronoi)] <- c("Bear1", "Bear2")
  
  # Get maximum distance of a sampled point (x_i, y_i) from the origin (0, 0):
  MaxDistanceFromOrigin <- max(sqrt((x ^ 2) + (y ^ 2)))
  
  # Generate plot radius from maximum distance from origin multiplied by some edge factor:
  PlotRadius <- EdgeFactor * MaxDistanceFromOrigin
  
  # Place heat values on a zero to one scale:
  RescaledHeat <- (heat - min(heat)) / max(heat - min(heat))
  
  # Set plot colour palette using viridis library (probably want to add options for this later):
  PlotColourPalette <- viridis::magma(n = NColours, begin = min(heat), end = max(heat))
  
  # Set plot colours by matching heat values to colour palette:
  PlotColours <- unlist(lapply(as.list(RescaledHeat), function(x) PlotColourPalette[max(which(x >= seq(from = 0, to = 1, length.out = NColours)))]))
  
  for(i in sort(c(which(DelaunayForVoronoi[, "bp1"]), which(DelaunayForVoronoi[, "bp2"])))) {
    
    # Get correct column (first or second):
    ColumnsToUse <- ifelse(DelaunayForVoronoi[i, "bp1"], "First", "Second")
    
    if(ColumnsToUse == "First") Coordinates <- DelaunayForVoronoi[i, c("x1", "y1")]
    if(ColumnsToUse == "Second") Coordinates <- DelaunayForVoronoi[i, c("x2", "y2")]
    
    x1 <- 0
    y1 <- 0
    x2 <- Coordinates[1]
    y2 <- Coordinates[2]
    
    Bearing <- (180 / pi) * atan2(as.numeric(y2 - y1), as.numeric(x2 - x1))
    
    Bearing <- Bearing %% 360
    
    x_new <- PlotRadius * cos(Bearing * (pi / 180))
    y_new <- PlotRadius * sin(Bearing * (pi / 180))
    
    if(ColumnsToUse == "First") DelaunayForVoronoi[i, c("x1", "y1")] <- c(x_new, y_new)
    if(ColumnsToUse == "Second") DelaunayForVoronoi[i, c("x2", "y2")] <- c(x_new, y_new)
    
  }
  
  
  
  
  LeftOutsideCirclePoints <- which(sqrt(apply(DelaunayForVoronoi[, c("x1", "y1")] ^ 2, 1, sum)) > PlotRadius)
  RightOutsideCirclePoints <- which(sqrt(apply(DelaunayForVoronoi[, c("x2", "y2")] ^ 2, 1, sum)) > PlotRadius)
  
  if(length(LeftOutsideCirclePoints) > 0) {
    
    NewX <- PlotRadius * cos(DelaunayForVoronoi[LeftOutsideCirclePoints, "Bear1"] * (pi / 180))
    NewY <- PlotRadius * sin(DelaunayForVoronoi[LeftOutsideCirclePoints, "Bear1"] * (pi / 180))
    
    DelaunayForVoronoi[LeftOutsideCirclePoints, c("x1", "y1")] <- cbind(NewX, NewY)
    
  }
  
  if(length(RightOutsideCirclePoints) > 0) {
    
    NewX <- PlotRadius * cos(DelaunayForVoronoi[RightOutsideCirclePoints, "Bear2"] * (pi / 180))
    NewY <- PlotRadius * sin(DelaunayForVoronoi[RightOutsideCirclePoints, "Bear2"] * (pi / 180))
    
    DelaunayForVoronoi[RightOutsideCirclePoints, c("x2", "y2")] <- cbind(NewX, NewY)
    
  }
  
  # Subfunction to get points that describe the convex hull surrounding each datapoint:
  GetPolygonPoints <- function(x, DelauneyMatrix, EdgeIndices, PlotRadius) {
    
    # Get points corresponding to the current index polygon:
    polygonpoints <- rbind(as.matrix(DelauneyMatrix[DelauneyMatrix[, "ind1"] == x, c("x1", "y1")]), as.matrix(DelauneyMatrix[DelauneyMatrix[, "ind2"] == x, c("x2", "y2")]))
    
    # If this is also an edge index (i.e., a polygon that has part of it's edge as the edge of the plot):
    if(sum(EdgeIndices == x) > 0) {
      
      # Set current index (edge point already in sample):
      CurrentIndex <- x
      
      # Find the rows of the Delauney matrix that correspond to the current index:
      IndexRows <- c(which(DelauneyMatrix[, "ind1"] == x), which(DelauneyMatrix[, "ind2"] == x))
      
      # Find adjacent index (will decide which additional point needs to be added:
      AdjacentIndex <- setdiff(c(DelauneyMatrix[IndexRows, ][which(DelauneyMatrix[IndexRows, "bp1"]), "ind1"], DelauneyMatrix[IndexRows, ][which(DelauneyMatrix[IndexRows, "bp2"]), "ind2"]), x)
      
      # Find the rows of the Delauney matrix that correspond to the adjacent (edge) index:
      AdjacentIndexRows <- c(which(DelauneyMatrix[, "ind1"] == AdjacentIndex), which(DelauneyMatrix[, "ind2"] == AdjacentIndex))
      
      # Build adjacent index block (just of rows corresponding to edge values):
      AdjacentIndexBlock <- matrix(c(as.vector(t(unname(DelauneyMatrix[AdjacentIndexRows, ][which(DelauneyMatrix[AdjacentIndexRows, "bp1"]), c("x1", "y1", "ind1"), drop = FALSE]))), as.vector(t(unname(DelauneyMatrix[AdjacentIndexRows, ][which(DelauneyMatrix[AdjacentIndexRows, "bp2"]), c("x2", "y2", "ind2"), drop = FALSE])))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("x", "y", "ind")))
      
      # Add additional (adjacent) edge point to polygon points
      polygonpoints <- rbind(polygonpoints, AdjacentIndexBlock[which(AdjacentIndexBlock[, "ind"] == AdjacentIndex), c("x", "y")])
      
      # Get bearing of current edge point:
      CurrentBearing <- c(DelauneyMatrix[intersect(which(DelauneyMatrix[, "ind1"] == CurrentIndex), which(DelauneyMatrix[, "bp1"])), "Bear1"], DelauneyMatrix[intersect(which(DelauneyMatrix[, "ind2"] == CurrentIndex), which(DelauneyMatrix[, "bp2"])), "Bear2"])
      
      # Get bearing of adjacent edge point:
      AdjacentBearing <- c(DelauneyMatrix[intersect(which(DelauneyMatrix[, "ind1"] == AdjacentIndex), which(DelauneyMatrix[, "bp1"])), "Bear1"], DelauneyMatrix[intersect(which(DelauneyMatrix[, "ind2"] == AdjacentIndex), which(DelauneyMatrix[, "bp2"])), "Bear2"])
      
      # If Adjacent bearing is true minimum bearing:
      if((AdjacentBearing - CurrentBearing) %% 360 >= (CurrentBearing - AdjacentBearing) %% 360) {
        
        # Set adjacent as minimum bearing:
        MinBearing <- AdjacentBearing
        
        # Set current as maximum bearing:
        MaxBearing <- CurrentBearing
        
      # If Current bearing is true minimum bearing:
      } else {
        
        # Set adjacent as maximum bearing:
        MaxBearing <- AdjacentBearing
        
        # Set current as minimum bearing:
        MinBearing <- CurrentBearing
        
      }
      
      # As long as there is space to add additional points to describe the curved edge of the circle:
      if(((MaxBearing - MinBearing) %% 360) > 1) {
        
        # Get interpolated bearings along edge of plot:
        InterpolatedBearings <- (ceiling(MinBearing) + c(0:((ifelse(floor(MaxBearing) == 0, 360, floor(MaxBearing)) - ceiling(MinBearing)) %% 360))) %% 360
        
        # Convert bearings to coordinates:
        InterpolatedCoordinates <- cbind(PlotRadius * cos(InterpolatedBearings * (pi / 180)), PlotRadius * sin(InterpolatedBearings * (pi / 180)))
        
        # Add interpolated points to polygon:
        polygonpoints <- rbind(polygonpoints, InterpolatedCoordinates)
        
      }
      
    }
    
    # Place polygon points in "order", i.e., using a convex hull:
    convexhullpoints <- chull(x = polygonpoints[, 1], y = polygonpoints[, 2])
    
    # Set order of polygon points using convex hull:
    polygonpoints <- polygonpoints[convexhullpoints, ]
    
    # Output polygon points ready for plotting:
    return(polygonpoints)
    
  }
  
  # Define deg indices (i.e., indices that share and edge with the edge of the plot):
  EdgeIndices <- sort(unname(unlist(c(DelaunayForVoronoi[, c("ind1", "ind2")]))[unlist(c(DelaunayForVoronoi[, c("bp1", "bp2")]))]))
  
  # Build polygon list (of indices, ready for plotting):
  PolygonsList <- lapply(as.list(1:length(x)), GetPolygonPoints, DelauneyMatrix = DelaunayForVoronoi, EdgeIndices = EdgeIndices, PlotRadius = PlotRadius)
  
  # Subset list as coordinates (creates room for adding additional data to each polygon, e.g., heat):
  PolygonsList <- lapply(PolygonsList, function(x) {x <- list(x); names(x) <- "PolygonCoordinates"; return(x)})
  
  # FOr each index:
  for(i in 1:N) {
    
    # Add plot colour to list:
    PolygonsList[[i]][[2]] <- PlotColours[i]
    
    # Add plot colour name to list:
    names(PolygonsList[[i]])[2] <- "PlotColour"
    
  }
  
  par(mar = rep(0, 4))
  
  # Create empty plot:
  plot(x = x, y = y, type = "n", xlim = c(-PlotRadius, PlotRadius), ylim = c(-PlotRadius, PlotRadius), asp = 1, axes = FALSE, xlab = "", ylab = "")
  
  # Add polygons coloured by heat values:
  Empty <- lapply(PolygonsList, function(x) polygon(x = x$PolygonCoordinates[, 1], y = x$PolygonCoordinates[, 2], col = x$PlotColour, border = 0))
  
}
