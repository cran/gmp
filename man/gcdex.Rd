\name{gcdex}
\alias{gcdex}

\title{Compute Bezoult coefficient}
\description{
  Compute g,s,t as
  \eqn{as + bt = g = gcd(a,b)}{as + bt = g = gcd(a,b)}.
  s and t are also known as Bezoult coefficients.
}
\usage{
gcdex(a, b)
}

\arguments{
  \item{a,b}{Either integer, numeric or string value
  (String value: ither starting with \code{0x} for hexadecimal, \code{0b} for
  binary or without prefix for decimal values.)
  Or an element of class bigz.}

}

\value{
  3 values:
  \item{g,s,t}{Elements of class bigz}

}

\references{The GNU MP Library, see \url{http://gmplib.org}}
\author{Antoine Lucas}


\seealso{\code{\link{gcd.bigz}}}
\examples{
gcdex(342,654)

}
\keyword{arith}

