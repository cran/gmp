\name{gmp.rand}
\alias{gmp.rand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a random number}
\description{
  Generate a uniformly distributed random number in the range 0 to
  \eqn{2^n -1}{2^n - 1}, inclusive.

  If seed is precise, a new seed is initialised.
}
\usage{
gmp.rand(n, seed = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ Integer}
  \item{seed}{Integer}
}


\value{
  A number.
}
\references{mpz\_urandomb from Gnu MP Library see http://swox.com/gmp}

\author{Antoine Lucas}



\examples{
# Integers are differents
gmp.rand(123)
gmp.rand(123)
gmp.rand(123)
# Integers are the same
gmp.rand(123,seed="234234234324323")
gmp.rand(123,seed="234234234324323")


}
\keyword{arith }% at least one, from doc/KEYWORDS

