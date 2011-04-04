\name{Wolff}
\alias{Wolff}
\title{Generate Random Samples from a Compound Potts Model by the Wolff Algorithm}
\description{
  Generate random samples from a compound Potts model using the Wolff Algorithm.
}
\usage{
  Wolff(n, nvertex, ncolor, neighbors, weights, beta)
}
\arguments{
  \item{n}{number of samples.}
  \item{nvertex}{number of vertices of a graph.}
  \item{ncolor}{number of colors each vertex can take.}
  \item{neighbors}{neighbors of a graph.}
  \item{weights}{weights between neighbors. One for each corresponding
	component in \code{neighbors}. The default values are 1s for all.}
  \item{beta}{the parameter inverse temperature of the Potts model.}
}
\value{
  A \code{nvertex} by \code{n} matrix with the kth column being the kth
  sample. 
}
\details{
  We use the Wolff algorithm to generate random samples from a compound
  Potts model. See \code{\link[PottsUtils]{rPotts1}} for more
  information on the compound Potts model.
}
\references{
 Ulli Wolff (1989)
 Collective Monte Carlo Updating for Spin Systems
 \emph{Physical Review Letters} \bold{vol. 62, no. 4}, 361-364
 
 Dai Feng (2008)
 Bayesian Hidden Markov Normal Mixture Models with Application to MRI
 Tissue Classification
 \emph{Ph. D. Dissertation, The University of Iowa} 
}
\seealso{
  \code{\link[PottsUtils]{SW}}, \code{\link[PottsUtils]{BlocksGibbs}}
}
\examples{
  #Example 1: Generate 100 samples from a Potts model with the
  #           neighborhood structure corresponding to a
  #           second-order Markov random field defined on a 
  #           3*3 2D graph. The number of colors is 2.
  #           beta=0.7. All weights are equal to 1.

  neighbors <- getNeighbors(mask=matrix(1, 3, 3), neiStruc=rep(2,4))
  Wolff(n=100, nvertex=9, ncolor=2, neighbors, beta=0.7)
}
\keyword{distribution}



