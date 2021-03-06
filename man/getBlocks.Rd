\name{getBlocks}
\alias{getBlocks}
\title{Get Blocks of a Graph}
\description{
  Obtain blocks of vertices of a 1D, 2D, or 3D graph, in order to use
  the conditional independence to speed up the simulation (checkerboard
  idea). 
}
\usage{
   getBlocks(mask, nblock)
}
\arguments{
  \item{mask}{a vector, matrix, or 3D array specifying vertices of a
    graph. Vertices of value 1 are within the graph and 0 are not.}
  \item{nblock}{a scalar specifying the number of blocks. For
    a 2D graph \code{nblock} could be either 2 or 4, and for a 3D
	graph \code{nblock} could be either 2 or 8.
    
  }
}
\value{
  A list with the number of components equal to \code{nblock}. Each
  component consists of vertices within the same block.
 }
\details{
  The vertices within each block are mutually independent given the
  vertices in other blocks. Some blocks could be empty.
}
\references{
    Darren J. Wilkinson
    Parallel Bayesian Computation
    \cite{Handbook of Parallel Computing and Statistics} 481-512
    Marcel Dekker/CRC Press  2005
}
\examples{
  #Example 1: split a line into 2 blocks
  getBlocks(mask=c(1,1,1,1,0,0,1,1,0), nblock=2)

  #Example 2: split a 4*4 2D graph into 4 blocks in order
  #           to use the checkerboard idea for a neighborhood structure
  #           corresponding to the second-order Markov random field.
  getBlocks(mask=matrix(1, nrow=4, ncol=4), nblock=4)

  #Example 3: split a 3*3*3 3D graph into 8 blocks
  #           in order to use the checkerboard idea for a neighborhood
  #           structure based on the 18 neighbors definition, where the
  #           neighbors of a vertex comprise its available
  #           adjacencies sharing the same edges or faces.
  mask <- array(1, dim=rep(3,3))
  getBlocks(mask, nblock=8)

}
\keyword{spatial}



