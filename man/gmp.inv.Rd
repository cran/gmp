\name{gmp.inv}
\alias{gmp.inv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute inverse of a modulo b.}
\description{
Compute inverse of a modulo b. If inverse does not exists, return 0.
}
\usage{
gmp.inv(a, b)
}

\arguments{
  \item{a}{Integer or string from an integer}
  \item{b}{Integer or string from an integer}
}

\references{Gnu MP Library see http://swox.com/gmp}


\author{Antoine Lucas}


\examples{
gmp.inv(343,464)

}
\keyword{arith}% at least one, from doc/KEYWORDS

