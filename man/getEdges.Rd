\name{getEdges}
\alias{getEdges}
\title{Get Edges of a Graph}
\description{
  Obtain edges of a 1D, 2D, or 3D graph based on the
  neighborhood structure.
}
\usage{
   getEdges(mask, neiStruc)
}
\arguments{
  \item{mask}{a vector, matrix, or 3D array specifying vertices of a
    graph. Vertices of value 1 are within the graph and 0 are not.}
  \item{neiStruc}{a scalar, vector of four components, or
	\eqn{3\times4}	matrix
    corresponding to 1D, 2D, or 3D graphs. It specifies the neighborhood
	structure. See \code{getNeighbors} for details.
    }
  
}
\value{
  A matrix of two columns with one edge per row.  The edges connecting
  vertices and their corresponding first neighbors are listed first, and
  then those corresponding to the second neighbors, and so on and so
  forth.  The order of neighbors is the same as in \code{getNeighbors}.
}
\details{
  There could be more than one way to define the same 3D neighborhood
  structure for a graph (see Example 4 for illustration). 
}
\references{
 Gerhard Winkler (1995)
 Image Analysis, Random Fields and Dynamic Monte Carlo Methods
 \emph{Springer-Verlag}

 Dai Feng (2008)
 Bayesian Hidden Markov Normal Mixture Models with Application to MRI
 Tissue Classification
 \emph{Ph. D. Dissertation, The University of Iowa} 

} 
\examples{
  #Example 1: get all edges of a 1D graph. 
  mask <- c(0,0,rep(1,4),0,1,1,0,0)
  getEdges(mask, neiStruc=2)

  #Example 2: get all edges of a 2D graph based on neighborhood structure
  #           corresponding to the first-order Markov random field.
  mask <- matrix(1 ,nrow=2, ncol=3)
  getEdges(mask, neiStruc=c(2,2,0,0))

  #Example 3: get all edges of a 2D graph based on neighborhood structure
  #           corresponding to the second-order Markov random field.
  mask <- matrix(1 ,nrow=3, ncol=3)
  getEdges(mask, neiStruc=c(2,2,2,2))

  #Example 4: get all edges of a 3D graph based on 6 neighbors structure
  #           where the neighbors of a vertex comprise its available
  #           N,S,E,W, upper and lower adjacencies. To achieve it, there
  #           are several ways, including the two below.
  mask <- array(1, dim=rep(3,3))
  n61 <- matrix(c(2,2,0,0,
                  0,2,0,0,
                  0,0,0,0), nrow=3, byrow=TRUE)
  n62 <- matrix(c(2,0,0,0,
                  0,2,0,0,
                  2,0,0,0), nrow=3, byrow=TRUE)
  e1 <- getEdges(mask, neiStruc=n61)
  e2 <- getEdges(mask, neiStruc=n62)
  e1 <- e1[order(e1[,1], e1[,2]),]
  e2 <- e2[order(e2[,1], e2[,2]),]
  all(e1==e2)

  #Example 5: get all edges of a 3D graph based on 18 neighbors structure
  #           where the neighbors of a vertex comprise its available
  #           adjacencies sharing the same edges or faces.
  #           To achieve it, there are several ways, including the one below.

  n18 <- matrix(c(2,2,2,2,
                  0,2,2,2,
                  0,0,2,2), nrow=3, byrow=TRUE)  
  mask <- array(1, dim=rep(3,3))
  getEdges(mask, neiStruc=n18)
  
  
}
\keyword{spatial}



