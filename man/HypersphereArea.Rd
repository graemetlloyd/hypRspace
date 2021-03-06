% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HypersphereArea.R
\name{HypersphereArea}
\alias{HypersphereArea}
\title{Get hypersphere area}
\usage{
HypersphereArea(Radius, NDimensions)
}
\arguments{
\item{Radius}{The radius of the sphere.}

\item{NDimensions}{The number of dimensions of the sphere.}
}
\value{
The surface area (or dimensional equivalent of area).
}
\description{
Calculates the dimensional equivalent of area for a hypersphere.
}
\details{
Given a specified radius and number of dimensions calculates the dimensional equivalent of area for a hypersphere.
}
\examples{

# Get surface area of the regular unit sphere:
HypersphereArea <- function(Radius = 1, NDimension = 3)

}
\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}
\keyword{hypersphere}
