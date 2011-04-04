\name{getPatches}
\alias{getPatches}
\title{Get Patches of a Graph}
\description{
  Obtain patches of a graph by Rem's algorithm.
}
\usage{
   getPatches(bonds, nvertex)
}
\arguments{
  \item{bonds}{a matrix of bonds in a graph, with one bond per row.}
  \item{nvertex}{number of vertices in a graph.}
}
\value{
  A list comprises all patches in a graph. Each component of the list
  consists of vertices within one patch.
}
\details{
  Given all bonds and the number of vertices in a graph,
  this function provides all patches. 
}
\references{
 Edsger W. Dijkstra (1976)
 A Discipline of Programming
 \emph{Englewood Cliffs, New Jersey : Prentice-Hall, Inc}
} 
\examples{
  #Example 1: Find patches of a 3*3 2D graph with 6 bonds.

  bonds <- matrix(c(1,2,2,5,5,6,3,6,5,8,7,8), ncol=2, byrow=TRUE)
  getPatches(bonds, 9)
}
\keyword{spatial}



