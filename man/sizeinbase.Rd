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
  \item{a}{Integer}
  \item{b}{base}
}

\value{
Integer: size (number of digits).
}

\references{Gnu MP Library see http://swox.com/gmp}
\author{Antoine Lucas}



\examples{
sizeinbase(342434,10)

}
\keyword{arith}

