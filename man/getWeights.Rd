\name{getWeights}
\alias{getWeights}
\title{Get All Weights of a Graph}
\description{
  Obtain weights of edges of a 1D, 2D, or 3D graph based on the
  neighborhood structure.
}
\usage{
   getWeights(mask, neiStruc, format=1)
}
\arguments{
  \item{mask}{a vector, matrix, or 3D array specifying vertices within a
    graph. Vertices of value 1 are within the graph and 0 are not.}
  \item{neiStruc}{a scalar, vector of four components,
	or \eqn{3\times4} matrix corresponding to 1D, 2D, or 3D graphs.
	It specifies the neighborhood structure.
	See \code{getNeighbors} for details.}
  \item{format}{If it is 1, then the output is a vector of weights,
	one for two vertices in the corresponding output from
	\code{getEdges}. If it is 2, then the output is a matrix, one for
	two vertices in the corresponding output from \code{getNeighbors}.
	The default value is 1.}
}
\value{
  A vector of weights, one component corresponding to an edge of a
  graph. Or a matrix of weights, one component corresponding to two
  vertices in neighbor.
}
\details{
  The weights are equal to the reciprocals of the distance between
  neighboring vertices.
}
\examples{
  #Example 1: get all weights of a 2D graph based on neighborhood structure
  #           corresponding to the first-order Markov random field.
  mask <- matrix(1 ,nrow=2, ncol=3)
  getWeights(mask, neiStruc=c(2,2,0,0))

  #Example 2: get all weights of a 2D graph based on neighborhood structure
  #           corresponding to the second-order Markov random field.
  #           Put the weights in a matrix form corresponding to
  #           neighbors of vertices.
  mask <- matrix(1 ,nrow=3, ncol=3)
  getWeights(mask, neiStruc=c(2,2,2,2), format=2)

  #Example 3: get all weights of a 3D graph based on 6 neighbors structure
  #           where the neighbors of a vertex comprise its available
  #           N,S,E,W, upper and lower adjacencies. 
  mask <- array(1, dim=rep(3,3))
  n61 <- matrix(c(2,2,0,0,
                  0,2,0,0,
                  0,0,0,0), nrow=3, byrow=TRUE)
  getWeights(mask, neiStruc=n61)
 }
\keyword{spatial}



