\name{modulus}
\alias{modulus}
\alias{modulus<-}
\alias{modulus.bigz}
\alias{modulus<-.bigz}

\title{Modulus}

\description{Access to modulus of a bigz}

\usage{
modulus(a,value)

}

\arguments{
  \item{a}{Object of class bigz}
  \item{value}{Object of class bigz}
}

\references{
   \item Gnu MP Library see \url{http://swox.com/gmp},
   \item Home page: \url{http://mulcyber.toulouse.inra.fr/projects/gmp/}
 }


\examples{
  x <- as.bigz(24)
  modulus(x)

  # x element of Z/31Z
  modulus(x) <- 31
  x+x
  # (6 %% 31)

  
}

\keyword{arith}
