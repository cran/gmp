\name{gcd.bigz}
\alias{gcd.bigz}
\alias{gcd}
\alias{gcd.default}
\alias{lcm.default}
\alias{lcm}
\alias{lcm.bigz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Great common divisor, Least common multiple}
\description{
  Great common divisor and least common multiple. 
}
\usage{
gcd.bigz(a, b)
lcm.bigz(a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a,b}{Either integer, numeric or string value
  (String value: ither starting with \code{0x} for hexadecimal, \code{0b} for
  binary or without prefix for decimal values.)
  Or an element of class bigz.}

}

\value{
  An element of class bigz
}
\references{Gnu MP Library see http://swox.com/gmp}
\author{Antoine Lucas}

\seealso{\code{\link{gcdex}}}

\examples{
gcd.bigz(210,342)
lcm.bigz(210,342)
}
\keyword{arith}% at least one, from doc/KEYWORDS

