\name{rPotts1}
\alias{rPotts1}
\title{Generate One Random Sample from a Potts Model}
\description{
  Generate a random sample from a Potts model with external field
  by Gibbs Sampling that takes advantage of conditional independence,
  or the partial decoupling method.
}
\usage{
  rPotts1(nvertex, ncolor, neighbors, blocks, edges=NULL, weights=1,
          spatialMat=NULL, beta, external, colors,
          algorithm=c("Gibbs", "PartialDecoupling"))
}
\arguments{
  \item{nvertex}{number of vertices in a graph.}
  \item{ncolor}{number of colors each vertex can take.}
  \item{neighbors}{all neighbors in a graph. It is not required 
	when using the partial decoupling method.}
\item{blocks}{the blocks of vertices in a graph.
  It is not required when using the partial decoupling method.}
  \item{edges}{all edges in a graph. The default value is \code{NULL}.
	It is not required 	when using Gibbs sampling.}
  \item{weights}{weights between neighbors or \eqn{\delta_{ij}}s
	in the partial decoupling method. When using Gibbs sampling,
	there is one for each corresponding
	neighbor in	\code{neighbors}. When using partial decoupling,
	there is one for each corresponding	edge in	\code{edges}. The default values are 1s for all.}
  \item{spatialMat}{a matrix that describes the relationship among
    vertices in neighbor. It is not required when using the partial decoupling method. The default value is \code{NULL} corresponding
	to the simple or compound Potts model.}
  \item{beta}{the parameter inverse temperature of the Potts model.}
  \item{external}{a matrix giving values of external field.  The number
	of rows equal to \code{nvertex} and number of columns equal to
	\code{ncolor}.}
  \item{colors}{the current colors of vertices.}
  \item{algorithm}{a character string specifying the algorithm used to
	generate samples. It must be either "Gibbs", or "PartialDecoupling",
	and may be abbreviated.	The default is "Gibbs".}
}
\value{
  The output is a vector with the kth component being the new color of
  vertex k.
}
\details{
  This function generates random samples from a 
  Potts model as follows:
  \deqn{
  p(z)=C(\beta)^{-1} \exp\{\sum_i \alpha_i(z_i) \beta \sum_{i \sim j} w_{ij} f(z_{i},z_{j})\}
  }

  where \eqn{C(\beta)} is a normalizing constant
  and \eqn{i \sim j}
  indicates neighboring vertices.  The parameter \eqn{\beta} is called the
  "inverse temperature", which determines the level of spatial homogeneity between
  neighboring vertices in the graph. We assume \eqn{\beta > 0}.
  The set \eqn{z_{1}, z_{2},
    \ldots,} comprises the indices to the colors of all vertices. Function
  \eqn{f(z_{i}, z_{j})} determines the relationship among vertices
  in neighbor. Parameter \eqn{w_{ij}} is the weight between vertex
  \eqn{i} and \eqn{j}. And \eqn{\sum_i \alpha_i(z_i)} is called the "external field".

  For the simple, the compound, and the simple repulsion Potts
  models,
  the external field is equal to 0.
  For the simple and the compound Potts model
  \eqn{f(z_{i}, z_{j}) = I(z_{i}=z_{j})}.
  Parameters \eqn{w_{ij}} are all equal for the simple Potts model but
  not so for the compound model.
  
  For the repulsion Potts model \eqn{f(z_{i}, z_{j}) = a_{1}}, if
  \eqn{z_{i} = z_{j}}; \eqn{f(z_{i}, z_{j}) = a_{2}}, if
  \eqn{|z_{i} - z_{j}|=1};  \eqn{f(z_{i}, z_{j}) = a_{3}} otherwise.
  
  
  The argument \code{spatialMat} is used to specify the
  relationship among vertices in neighbor. The default value is
  \code{NULL}
  corresponding to the simple or the compound Potts 
  model.  The component at the \eqn{i}th row and \eqn{j}th column 
  defining the relationship when the color of a vertex is \eqn{i} and the
  color of its neighbors is \eqn{j}.
  Besides the default setup, for the simple and the compound Potts models
  \code{spatailMat} could be an identity matrix also. For
  the repulsion Potts model, it is
  \deqn{\left(\begin{array}{ccccc}
	a_1 & a_2 & a_3 & \ldots &  a_3 \\
	a_2 & a_1 & a_2 & \ldots &  a_3 \\
	\vdots & \vdots & \vdots & \ddots  & \vdots \\  
	a_3 & a_3 & a_3 & \ldots &  a_1
	\end{array}\right)}
  Other relationships among neighboring vertices can be specified through
  it as well.

  The argument \code{external} gives the values of external field.
   
  Gibbs sampling can be used to generate samples from all
  kinds of Potts models. We use the method that takes advantage of
  conditional independence to speed up the simulation. See
  \code{\link[PottsUtils]{BlocksGibbs}} for details.

  The partial decoupling method could be used to generate samples
  from the simple Potts model plus the external field.
  The \eqn{\delta_{ij}}s are specified through the argument \code{weights}.
  
}
\references{
 Dai Feng (2008)
 Bayesian Hidden Markov Normal Mixture Models with Application to MRI
 Tissue Classification
 \emph{Ph. D. Dissertation, The University of Iowa}

 David M. Higdon (1998)
 Auxiliary variable methods for Markov Chain Monte Carlo with
 applications
 \emph{Journal of the American Statistical Association}
 \bold{vol. 93} 585-595
} 
\seealso{
  \code{\link[PottsUtils]{BlocksGibbs}}, \code{\link[PottsUtils]{Wolff}}
  \code{\link[PottsUtils]{SW}}
}
\examples{
  \dontrun{
  neighbors <- getNeighbors(matrix(1, 16, 16), c(2,2,0,0))
  blocks <- getBlocks(matrix(1, 16, 16), 2)
  spatialMat <- matrix(c(2, 0, -1, 0, 2, 0, -1, 0, 2), ncol=3)
  mu <- c(22, 70 ,102)
  sigma <- c(17, 16, 19)
  count <- c(40, 140, 76)
  y <- unlist(lapply(1:3, function(i) rnorm(count[i], mu[i], sigma[i])))
  external <- do.call(cbind,
                      lapply(1:3, function(i) dnorm(y, mu[i],sigma[i])))
  current.colors <- rep(1:3, count)
  rPotts1(nvertex=16^2, ncolor=3, neighbors=neighbors, blocks=blocks,  
          spatialMat=spatialMat, beta=0.3, external=external,
          colors=current.colors, algorithm="G")
  edges <- getEdges(matrix(1, 16, 16), c(2,2,0,0))
  rPotts1(nvertex=16^2, ncolor=3, edges=edges, beta=0.3,
          external=external, colors=current.colors, algorithm="P")
 }
}
\keyword{distribution}
