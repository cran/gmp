\name{gcd.bigz}
\alias{gcd.bigz}
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
  \item{a,b}{Object of class bigz}xs

}

\value{
  A string of the integer value.
}
\references{
   \item Gnu MP Library see \url{http://swox.com/gmp},
   \item Home page: \url{http://mulcyber.toulouse.inra.fr/projects/gmp/}
 }


\author{Antoine Lucas}

\seealso{\code{\link{gcdex}}}

\examples{
gcd.bigz(210,342)
lcm.bigz(210,342)
}
\keyword{arith}% at least one, from doc/KEYWORDS

