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
  \item{n}{Integer}
}
\details{
  This function uses probabilistic algorithm to identify
  primes. For practical purposes, it's adequat, the chance
  of a composit passing will be extremely small.
}
\value{
  A (probably) prime  number
}
\references{Gnu MP Library see http://swox.com/gmp}

\author{Antoine Lucas}


\seealso{\code{\link{isprime}}}
\examples{
nextprime(14)
}
\keyword{arith}% at least one, from doc/KEYWORDS


