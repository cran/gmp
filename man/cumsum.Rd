\name{cumsum}
\alias{cumsum.bigz}
\alias{cumsum.bigq}
\alias{sum.bigz}
\alias{sum.bigq}
\alias{prod}
\alias{prod.bigz}
\alias{prod.bigq}
\alias{prod.default}

\title{Cummulative, product}

\usage{
cumsum.bigz(...)
cumsum.bigq(...)
sum.bigz(...)
sum.bigq(...)
prod(...)
}

\description{
  Theses functions overload sum, cumsum and prod function for big rationals and
  big integers
}

\arguments{
  \item{...}{an element of class bigz or bigq.}
  }
  

\value{
  return an element of class bigz or bigq.
}

\author{Antoine Lucas}

\examples{
  x = as.bigz(1:12)
  cumsum(x)
  prod(x)
  sum(x)

  x = as.bigq(1:12)
  cumsum(x)
  prod(x)
  sum(x)
}

\keyword{arith}


\seealso{\code{\link[base]{apply}}}

