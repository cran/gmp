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
  \item{a}{Integer or string from an integer}
  \item{b}{Integer or string from an integer}
}

\value{
  \item{g}{Integer}
  \item{s}{Integer}
  \item{t}{Integer}

}

\references{Gnu MP Library see http://swox.com/gmp}
\author{Antoine Lucas}


\seealso{\code{\link{gcd.bigz}}}
\examples{
gcdex(342,654)

}
\keyword{arith}

