\name{Extract}

\alias{[.bigz}
\alias{[<-.bigz}
\alias{[[.bigz}
\alias{[[<-.bigz}
\alias{c.bigz}
\alias{rep.bigz}
\alias{length.bigz}
\alias{length<-.bigz}
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
\method{c}{bigz}(..., recursive = FALSE)
\method{c}{bigq}(..., recursive = FALSE)
\method{rep}{bigz}(x,times, ...)
\method{rep}{bigq}(x,times, ...)
}

\arguments{
  \item{x}{\R object of class \code{"bigz"} or \code{"bigq"}, respectively.}
  \item{...}{further arguments, notably for \code{c()}.}
  \item{times}{integer}
  \item{recursive}{unused here}
}

\note{
  Unlike standard matrices, \code{x[i]} and \code{x[i,]} do the same.
  %% MM: hmm, just because of lazyness ??
}

\examples{
  a <- as.bigz(123)
  ## indexing "outside" --> extends the vectors (filling with NA)
  a[2] <- a[1]
  a[4] <- -4

  ## create a vector of 3 a
  c(a,a,a)

  ## repeate a 5 times
  rep(a,5)

  ## with matrix
  m <- matrix.bigz(1:6,3)

  ## these do the same:
  m[1,]
  m[1]
  m[-c(2,3),]
  m[-c(2,3)]
  m[c(TRUE,FALSE,FALSE)]
\dontshow{
  stopifnot(identical(c(a,a), rep(a,2)),
            identical(rep(a,0), a[FALSE]), length(a[FALSE]) == 0,
            identical(m[1,], m[-(2:3),]),
            identical(m[1,], m[c(TRUE,FALSE,FALSE),])
  )
}
  ##_modification on matrix
  m[2,-1] <- 11
}
\keyword{arith}

