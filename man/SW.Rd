\name{SW}
\alias{SW}
\title{Generate Random Samples from a Compound Potts Model by the Swendsen-Wang Algorithm}
\description{
  Generate random samples from a compound Potts model using the Swendsen-Wang algorithm.
}
\usage{
  SW(n, nvertex, ncolor, edges, weights, beta)
}
\arguments{
  \item{n}{number of samples.}
  \item{nvertex}{number of vertices of a graph.}
  \item{ncolor}{number of colors each vertex can take.}
  \item{edges}{edges of a graph.}
  \item{weights}{weights of edges. One for each corresponding component in
	\code{edges}. The default values are 1s for all.}
  \item{beta}{the parameter inverse temperature of the Potts model.}
}
\value{
  The output is a \code{nvertex} by \code{n} matrix with
  the kth column being the kth sample. 
}
\details{
  We use the Swendsen-Wang algorithm to generate random samples from a
  compound Potts model. See \code{\link[PottsUtils]{rPotts1}} for more
  information on the compound Potts model.
 }
\references{
 Robert H. Swendsen and Jian-Sheng Wang (1987)
 Nonuniversal Critical Dynamics in Monte Carlo Simulations
 \emph{Physical Review Letters} \bold{vol. 58, no. 2}, 86-88

 Dai Feng (2008)
 Bayesian Hidden Markov Normal Mixture Models with Application to MRI
 Tissue Classification
 \emph{Ph. D. Dissertation, The University of Iowa} 
} 
\seealso{
  \code{\link[PottsUtils]{Wolff}}, \code{\link[PottsUtils]{BlocksGibbs}}
}
\examples{
  #Example 1: Generate 100 samples from a Potts model with the
  #           neighborhood structure corresponding to a
  #           second-order Markov random field defined on a 
  #           3*3 2D graph. The number of colors is 2.
  #           beta=0.1. All weights are equal to 1.
 
  edges <- getEdges(mask=matrix(1, 3, 3), neiStruc=rep(2,4))
  SW(n=100, nvertex=9, ncolor=2, edges, beta=0.1)
}
\keyword{distribution}
