\name{getNC}
\alias{getNC}
\title{Calculate the Normalizing Constant of a Simple Potts Model}
\description{
  Use the thermodynamic integration approach to calculate the normalizing constant of a Simple Potts Model.
}
\usage{
  getNC(beta, subbetas, nvertex, ncolor,
        edges, neighbors=NULL, blocks=NULL, 
        algorithm=c("SwendsenWang", "Gibbs", "Wolff"), n, burn)
}
\arguments{
  \item{beta}{the inverse temperature parameter of the Potts model.}
  \item{subbetas}{vector of betas used for the integration.}
  \item{nvertex}{number of vertices in a graph.}
  \item{ncolor}{number of colors each vertex can take.}
  \item{edges}{all edges in a graph.}
  \item{neighbors}{all neighbors in a graph. The default is
    \code{NULL}. If the sampling algorithm is "BlocksGibbs" or "Wolff",
    then this has to be specified.}
  \item{blocks}{the blocks of vertices of a graph. The default is
    \code{NULL}.  If the sampling algorithm is "BlocksGibbs", 
    then this has to be specified.}
  \item{algorithm}{a character string specifying the algorithm used to
	generate samples. It must be one of "SwendsenWang", "Gibbs", or
	"Wolff" and may be abbreviated. The default is "SwendsenWang".}
  \item{n}{number of iterations.}
  \item{burn}{number of burn-in.}

}
\value{
  The corresponding normalizing constant.
}
\details{
  Use the thermodynamic integration approach to calculate the normalizing
  constant from a simple Potts model. See
  \code{\link[PottsUtils]{rPotts1}}
  for more information on the simple Potts model. 
  

  By the thermodynamic integration method,
  \deqn{
   \log{C(\beta)} =  N\log{k} + \int_{0}^{\beta}E(U({\bf z})|\beta^{'}, k)d\beta^{'}
  }
  where N is the total number of vertices (\code{nvertex}), k is the
  number of colors (\code{ncolor}), and
  \eqn{U({\bf z}) = \sum_{i \sim j}\textrm{I}(z_{i}=z_{j})}.
  Calculate \eqn{E(U({\bf z})} for \code{subbetas} based on samples,
  and then compute the integral by numerical integration.
}
\references{
 Peter J. Green and Sylvia Richardson  (2002)
 Hidden Markov Models and Disease Mapping 
 \emph{Journal of the American Statistical Association} \bold{vol. 97, no. 460}, 1055-1070
} 
\seealso{
  \code{\link[PottsUtils]{BlocksGibbs}}, \code{\link[PottsUtils]{SW}},
  \code{\link[PottsUtils]{Wolff}}
}
\examples{
  #Example 1: Calculate the normalizing constant of a simple Potts model
  #           with the neighborhood structure corresponding to a
  #           first-order Markov random field defined on a
  #           \eqn{3 \times3} 2D graph. The number of colors is 2 and
  #           \eqn{\beta=2}.
  #           Use 11 subbetas evenly distributed between 0 and 2.
  #           The sampling algorithm is Swendsen-Wang with 10000
  #           iterations and 1000 burn-in. 
 
  edges <- getEdges(mask=matrix(1,3,3), neiStruc=c(2,2,0,0))
  getNC(beta=2, subbetas=seq(0,2,by=0.2), nvertex=3*3, ncolor=2,
        edges, algorithm="S", n=10000, burn=1000)

}
\keyword{distribution}
