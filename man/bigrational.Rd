\name{bigq}
\alias{bigq}
\alias{as.bigq}
\alias{as.character.bigq}
\alias{as.double.bigq}
\alias{print.bigq}
\alias{is.na.bigq}
\alias{as.bigz.bigq}
\alias{denominator}
\alias{denominator<-}
\alias{numerator}
\alias{numerator<-}
\title{Large sized rationals}
\description{
  Type class supporting arithmetic operations on very large rationals.
}
\usage{
as.bigq(n, d = 1)
as.character.bigq(a)
as.double.bigq(x,...)
as.bigz.bigq(a,mod=NA)
is.na.bigq(a)
print.bigq(x,...)
denominator(x)
numerator(x)
}

\arguments{

\item{a,n,d}{Either integer, numeric or string value
  (String value: ither starting with \code{0x} for hexadecimal, \code{0b} for
  binary or without prefix for decimal values.
  Any format error results in \code{0})}
\item{mod}{Optional modulus to convert into biginteger}
\item{x}{Numeric value}
\item{...}{Additional parameters}

}
\value{
  A bigq class representing the parameter value.
}
\details{

}

\references{http://mulcyber.toulouse.inra.fr/projects/gmp/}


\author{Antoine Lucas}
%\note{ ~~further notes~~ }



\examples{
x <- as.bigq(21,6)
x
# 7 / 2
# Wow ! result is simplify.

y <- as.bigq(5,3)

# addition work !
x + y

# You can even try multiplication, division...
x * y / 13 

# convert to string, double
as.character(x)


}

\keyword{arith}
%\keyword{bigz}

