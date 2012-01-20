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
\S3method{as.character}{bigq}(x, b=10,...)
\S3method{as.double}{bigq}(x,...)
as.bigz.bigq(a, mod=NA)
\S3method{is.na}{bigq}(x)
\S3method{print}{bigq}(x, quote=FALSE, ...)
denominator(x)
numerator(x)
}

\arguments{
  \item{n,d}{either integer, numeric or string value
    (String value: either starting with \code{0x} for hexadecimal, \code{0b} for
    binary or without prefix for decimal values.
    Any format error results in \code{0}).
    \code{n} stands for numerator, \code{d} for denominator}
  \item{a}{an element of class \code{"bigq"}}
  \item{mod}{optional modulus to convert into biginteger}
  \item{x}{numeric value}
  \item{b}{base: from 2 to 36}
  \item{...}{additional arguments passed to methods}
  \item{quote}{(for printing:) logical indicating if the numbers
    should be quoted (as characters are);  the default used to be
    \code{TRUE} (implicitly) till 2011.}
}
\value{
  An \R object of (S3) class \code{"bigq"} representing the parameter value.
}
\details{
  \code{as.bigz.bigq()} returns the smallest integers not less than the
  corresponding rationals bigq.
}

\references{\url{http://mulcyber.toulouse.inra.fr/projects/gmp/}}

\author{Antoine Lucas}

\examples{
x <- as.bigq(21,6)
x
# 7 / 2
# Wow ! result is simplified.

y <- as.bigq(5,3)

# addition works !
x + y

# You can even try multiplication, division...
x * y / 13

# convert to string, double
as.character(x)
as.double(x)
}

\keyword{arith}
%\keyword{bigz}

