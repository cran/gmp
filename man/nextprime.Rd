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
\references{The GNU MP Library, see \url{http://gmplib.org}}

\author{Antoine Lucas}


\seealso{\code{\link{isprime}}}
\examples{
nextprime(14)
}
\keyword{arith}% at least one, from doc/KEYWORDS


