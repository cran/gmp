\name{nextprime}
\alias{nextprime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Next prime number}
\description{
Return next prime number geater than n
}
\usage{
nextprime(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{An integer, numeric, string or bigz}
}
\details{
  This function uses probabilistic algorithm to identify
  primes. For practical purposes, it's adequat, the chance
  of a composit passing will be extremely small.
}
\value{
  A (probably) prime  number
}

\references{
   \item Gnu MP Library see \url{http://swox.com/gmp},
   \item Home page: \url{http://mulcyber.toulouse.inra.fr/projects/gmp/}
 }

\author{Antoine Lucas}


\seealso{\code{\link{isprime}}}
\examples{
nextprime(14)
}
\keyword{arith}% at least one, from doc/KEYWORDS


