\name{getConfs}
\alias{getConfs}
\title{Generate Configurations of a Graph}
\description{
  Using recursive method to generate all possible configurations of a graph.
}
\usage{
   getConfs(nvertex, ncolor)
}
\arguments{
  \item{nvertex}{number of vertices in a graph.}
  \item{ncolor}{number of colors each vertex can take.}
}
\value{
  A matrix of all possible configurations. Each column corresponds
  to one configuration.
}
\details{
  Suppose there are n vertices and each can take values from
  \eqn{1,2, \ldots, ncolor}.
  This function generates all possible configurations. 
  For example, if there are two vertices and each can be either 1 or 2,
  then the possible configurations are (1,1), (1,2), (2,1) and (2,2).
}
\examples{
  #Example 1: There are two vertices and each is either of
  #           color 1 or 2.
  getConfs(2,2)
}
\keyword{math}



