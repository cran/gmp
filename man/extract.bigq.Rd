\name{Extract}

\alias{[.bigq}
\alias{[<-.bigq}
\alias{[[.bigq}
\alias{[[<-.bigq}
\alias{c.bigq}
\alias{rep.bigq}
\alias{length.bigq}
\alias{length<-.bigq}
\title{Extract or Replace Parts of an Object}
\description{
  Operators acting on vectors, arrays and lists to extract or replace subsets.
}
\usage{
c.bigq(...,recursive = FALSE)
rep.bigq(x,times,...)
}

\arguments{
  \item{x}{Object or vector of class bigz}
%  \item{i}{Integer, index of the vector}
  \item{...}{Additional parameters}
  \item{times}{Integer}
  \item{recursive}{Unused}
}

\note{
  Unlike standard matrix, operator x[i] & x[i,] do the same.

  }


\examples{
  a <- as.bigq(123)
  a[2] <- a[1]

  # create a vector of 3 a
  c(a,a,a)

  # repeate a 5 times
  rep(a,5)

  }

\keyword{arith}

