\name{==.bigz}

\alias{==.bigz}
\alias{!=.bigz}
\alias{<.bigz}
\alias{<=.bigz}
\alias{>.bigz}
\alias{>=.bigz}
\title{Relational Operators}
\description{
Binary operators which allow the comparison of values in atomic vectors.
}
\usage{
a < b
a > b
a <= b
a >= b
a == b
a != b
}

\arguments{
  \item{a,b}{Object or vector of class bigz}
}

\references{
   \item Gnu MP Library see \url{http://swox.com/gmp},
   \item Home page: \url{http://mulcyber.toulouse.inra.fr/projects/gmp/}
 }


\examples{
  x <- as.bigz(8000)
  x ^ 300 < 2 ^x
  }

\keyword{arith}

