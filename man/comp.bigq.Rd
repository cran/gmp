\name{Bigq}

\alias{==.bigq}
\alias{!=.bigq}
\alias{<.bigq}
\alias{<=.bigq}
\alias{>.bigq}
\alias{>=.bigq}
\alias{sign.bigq}
\title{Relational Operators}
\description{
  Binary operators which allow the comparison of values in atomic vectors.
}
\usage{
sign.bigq(a)
}

\arguments{
  \item{a}{Object or vector of class bigq}
}

\examples{
  x <- as.bigq(8000,21)
  x  < 2 * x
  }

\keyword{arith}

