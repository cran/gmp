\name{sizeinbase}
\alias{sizeinbase}

\title{Compute size of a biginteger in a base}
\description{
  Return size (number of digits) written in base b. 
}
\usage{
sizeinbase(a, b=10)
}

\arguments{
  \item{a}{An integer, numeric, string or bigz}
  \item{b}{base}
}

\value{
Integer: size (number of digits).
}


\references{
   \item Gnu MP Library see \url{http://swox.com/gmp},
   \item Home page: \url{http://mulcyber.toulouse.inra.fr/projects/gmp/}
 }


\author{Antoine Lucas}



\examples{
sizeinbase(342434,10)

}
\keyword{arith}

