\name{Relational Operator}

\alias{==.bigz}
\alias{!=.bigz}
\alias{<.bigz}
\alias{<=.bigz}
\alias{>.bigz}
\alias{>=.bigz}
\alias{sign.bigz}
\title{Relational Operators}
\description{
Binary operators which allow the comparison of values in atomic vectors.
}
\usage{
sign.bigz(a)
}

\arguments{
  \item{a}{Object or vector of class bigz}
}

\examples{
  x <- as.bigz(8000)
  x ^ 300 < 2 ^x
  }

\keyword{arith}

