\name{rMultimon}
\alias{rMultinom}
\title{Generate Random Samples from Different Multinomial
  Distributions with the Same Number of Classes
}
\description{
  Generate random samples from multinomial distributions with different
  classes probabilities.
}
\usage{
   rMultinom(p)
}
\arguments{
  \item{p}{matrix with each row specifying the probabilities for k classes
           of each sample.}
}
\value{
  Random samples from multinomial distributions.
}
\details{
  This function victories the generation of random samples from
  different multinomial distributions by the inversion of CDF method.
}
\seealso{
  \code{\link[stats]{rmultinom}}
}
\examples{
  #Example 1: Generate 100000 random samples from multinomial distribution
  #           with 3 classes all having the same probabilities.
  x <- rMultinom(p=matrix(1/3, nrow=100000, ncol=3))
  hist(x)
}
\keyword{distribution}



