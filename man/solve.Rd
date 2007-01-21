\name{solve.bigz}
\alias{solve.bigz}
\alias{solve.bigq}

\title{Solve a system of equation}
\description{
  This generic function solves the equation 'a \%*\% x = b' for 'x',
  where 'b' can be either a vector or a matrix.

  If a and b are rational, return is a rational matrix.

  If a and b are big integers (of class bigz) solution is in Z/nZ if
  there is a common modulus, or a rational matrix if not.
}
\usage{
solve.bigz(a, b, ...)
solve.bigq(a, b, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a,b}{A element of class bigz or bigq}
  \item{\dots}{Unused}
}
\details{
  It uses the Gauss and trucmuch algo \dots to be detailled.
}
\value{
  If a and b are rational, return is a rational matrix.

  If a and b are big integers (of class bigz) solution is in Z/nZ if
  there is a common modulus, of a rational matrix if not.

}

\author{Antoine Lucas}


\seealso{\code{\link[base]{solve}}}

\examples{
 x = matrix(1:4,2,2)
 ## standard solve
 solve(x)

 q = as.bigq(x)
 ## solve with rational
 solve(q)

 z = as.bigz(x)
 modulus(z) <- 7
 ## solve in Z/7Z
 solve(z)

 b = c(1,3)

 solve(q,b)

 solve(z,b)
}
\keyword{arith}% at least one, from doc/KEYWORDS

