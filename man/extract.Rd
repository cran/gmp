\name{[.bigz}

\alias{[.bigz}
\alias{[<-.bigz}
\alias{c.bigz}
\alias{rep.bigz}
\alias{length.bigz}
\alias{length<-.bigz}
\title{Extract or Replace Parts of an Object}
\description{
Operators acting on vectors, arrays and lists to extract or replace subsets.
}
\usage{
a[i]
c.bigz(...,recursive = FALSE)
rep.bigz(x,times,...)
}

\arguments{
  \item{a,x}{Object or vector of class bigz}
  \item{i}{Integer, index of the vector}
  \item{...}{Additional parameters}
  \item{times}{Integer}
  \item{recursive}{Unused}
}

\references{
   \item Gnu MP Library see \url{http://swox.com/gmp},
   \item Home page: \url{http://mulcyber.toulouse.inra.fr/projects/gmp/}
 }

\examples{
  a <- as.bigz(123)
  a[2] <- a[1]

  # create a vector of 3 a
  c(a,a,a)

  # repeate a 5 times
  rep(a,5)

  }

\keyword{arith}

