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
  \item{a}{An integer, numeric, string or bigz}
  \item{b}{An integer, numeric, string or bigz}
}

\value{
  \item{g}{Integer}
  \item{s}{Integer}
  \item{t}{Integer}

}

\references{
   \item Gnu MP Library see \url{http://swox.com/gmp},
   \item Home page: \url{http://mulcyber.toulouse.inra.fr/projects/gmp/}
 }

\author{Antoine Lucas}


\seealso{\code{\link{gcd.bigz}}}
\examples{
gcdex(342,654)

}
\keyword{arith}

