\name{BlocksGibbs}
\alias{BlocksGibbs}
\title{Generate Random Samples from a Potts Model Using the Checkerboard Idea
 }
\description{
  Generate random samples from a Potts model by Gibbs Sampling
  that takes advantage of conditional independence.
}
\usage{
  BlocksGibbs(n, nvertex, ncolor, neighbors, blocks,
              weights=1, spatialMat=NULL, beta)
}
\arguments{
  \item{n}{number of samples.}
  \item{nvertex}{number of vertices in a graph.}
  \item{ncolor}{number of colors each vertex can take.}
  \item{neighbors}{a matrix of all neighbors in a graph, one row per vertex.}
  \item{blocks}{a list of blocks of vertices in a graph.}
  \item{weights}{weights between neighbors. One for each corresponding
	neighbor in \code{neighbors}. The default values are 1s for all.}
  \item{spatialMat}{the matrix that describes the relationship among
    vertices in neighbor. The default value is \code{NULL} corresponding
	to the simple or compound Potts model.}
  \item{beta}{the parameter inverse temperature of the Potts model.}
 }
\value{
  The output is a \code{nvertex}*\code{n} matrix with
  the kth column being the kth sample. 
}
\details{
  We use the Gibbs algorithm that takes advantage of  
  conditional independence to speed up the generation of random samples from a 
  Potts model. The idea is that if we can divide variables that need to
  be updated into different blocks and given the variables in other
  blocks, all the variables within the same block are conditionally
  independent, then we can update all blocks iteratively with the variables
  within the same block being updated simultaneously.

  The \code{spatialMat} is the argument used to specify the
  relationship among vertices in neighbor. See
  \code{\link[PottsUtils]{rPotts1}} for more information on the Potts
  model and \code{spatialMat}.
}
\references{
 Dai Feng (2008)
 Bayesian Hidden Markov Normal Mixture Models with Application to MRI
 Tissue Classification
 \emph{Ph. D. Dissertation, The University of Iowa} 
} 
\seealso{
  \code{\link[PottsUtils]{Wolff}}, \code{\link[PottsUtils]{SW}}
}
\examples{
  #Example 1: Generate 100 samples from a repulsion Potts model with the
  #           neighborhood structure corresponding to a first-order
  #           Markov random field defined on a \eqn{3 \times 3} 2D graph.
  #           The number of colors is 3 and \eqn{\beta=0.1,a_1=2,a_2=1,a_3=0}.
  #           All weights are equal to 1.
 
  neighbors <- getNeighbors(mask=matrix(1, 3, 3), neiStruc=c(2,2,0,0))
  blocks <- getBlocks(mask=matrix(1, 3, 3), nblock=2)
  spatialMat <- matrix(c(2,1,0, 1,2,1,0,1,2), ncol=3)
  BlocksGibbs(n=100, nvertex=9, ncolor=3, neighbors=neighbors, blocks=blocks,
              spatialMat=spatialMat, beta=0.1)
}
\keyword{distribution}
