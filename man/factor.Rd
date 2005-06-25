\name{factorization}
\alias{factorize}

\title{Factorize a number}

\description{Give all primes numbers to factor the number}

\usage{factorize(n)}

\arguments{
  \item{n}{Either integer, numeric or string value
  (String value: ither starting with \code{0x} for hexadecimal, \code{0b} for
  binary or without prefix for decimal values.)
  Or an element of class bigz.}
}
\details{
  The factorization function uses the Pollard Rho algorithm.
  }

\value{
  Vector of class bigz.
}

\references{Gnu MP Library see http://swox.com/gmp}
\author{Antoine Lucas}



\examples{
factorize(34455342)
}

\keyword{arith}

