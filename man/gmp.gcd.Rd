\name{gcd}
\alias{gcd}
\alias{lcm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Great common divisor, Least common multiple}
\description{
  Great common divisor and least common multiple. In french: pgcd and
  ppcm.
}
\usage{
gcd(a, b)
lcm(a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{Integer or string from an integer}
  \item{b}{Integer or string from an integer}

}

\value{
  A string of the integer value.
}
\references{Gnu MP Library see http://swox.com/gmp}
\author{Antoine Lucas}

\seealso{\code{\link{gcdex}}}

\examples{
gcd(210,342)
lcm(210,342)
}
\keyword{arith}% at least one, from doc/KEYWORDS

