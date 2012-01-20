\name{gcd.bigz}
\title{Greatest Common Divisor, Least Common Multiple}
\alias{gcd.bigz}
\alias{gcd}
\alias{gcd.default}
\alias{lcm.default}
\alias{lcm}
\alias{lcm.bigz}
\description{
  Compute the greatest common divisor (GCD) and least common multiple
  (LCD) of two (big) integers.
}
\usage{
\method{gcd}{bigz}(a, b)
lcm.bigz(a, b)
}
\arguments{
  \item{a,b}{Either integer, numeric, \code{\link{bigz}} or a string
    value; if a string, either starting with \code{0x} for hexadecimal,
    \code{0b} for binary or without prefix for decimal values.}
}
\value{
  An element of class bigz
}
\references{The GNU MP Library, see \url{http://gmplib.org}}
\author{Antoine Lucas}

\seealso{\code{\link{gcdex}}}

\examples{
gcd.bigz(210,342) # or also
lcm.bigz(210,342)
a <- 210 ; b <- 342
stopifnot(gcd.bigz(a,b) * lcm.bigz(a,b) == a * b)
}
\keyword{arith}

