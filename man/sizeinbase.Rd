\name{sizeinbase}
\alias{sizeinbase}

\title{Compute size of a bigz in a base}
\description{
  Return size (number of digits) written in base b.
}
\usage{
sizeinbase(a, b=10)
}

\arguments{
  \item{a}{big integer, i.e. \code{"\link{bigz}"}}
  \item{b}{base}
}

\value{
  integer of the same length as \code{a}: the size, i.e. number of
  digits, of each \code{a[i]}.
}
\references{The GNU MP Library, see \url{http://gmplib.org}}
\author{Antoine Lucas}

\examples{
sizeinbase(342434, 10)# 6 obviously

Iv <- as.bigz(2:7)^500
sizeinbase(Iv)
stopifnot(sizeinbase(Iv)       == nchar(as.character(Iv)),
          sizeinbase(Iv, b=16) == nchar(as.character(Iv, b=16)))
}
\keyword{arith}

