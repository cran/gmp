\name{powm}
\alias{powm}
\alias{pow}

\title{Exponential function}
\description{
  powm compute \eqn{a ^ exp [modulo]}{a ^ exp [modulo]}
  pow compute \eqn{a ^ exp}{a ^ exp}
}
\usage{
powm(a, exp, modulo)
pow(a, exp)
}

\arguments{
  \item{a,exp,modulo}{Integer}
}

\value{
  Integer
}
\references{Gnu MP Library see http://swox.com/gmp}


\author{Antoine Lucas}


\examples{
powm(23,23,454545)
pow(1234,5)
1234^5
}
\keyword{arith}

