\name{bigq}
\alias{bigq}
\alias{bigq-class}
\alias{as.bigq}
\alias{is.bigq}
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
  Class \code{"bigq"} encodes rationals encoded as ratios of arbitrary
  large integers (via GMP).
  A simple S3 class (internally a \code{\link{raw}} vector), it has been
  registered as formal (S4) class (via \code{\link{setOldClass}}), too.

}
\usage{
as.bigq(n, d = 1)
\S3method{as.character}{bigq}(x, b=10,...)
\S3method{as.double}{bigq}(x,...)
as.bigz.bigq(a, mod=NA)
is.bigq(x)
\S3method{is.na}{bigq}(x)
\S3method{print}{bigq}(x, quote=FALSE, initLine = TRUE, ...)
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
  \item{x}{a \dQuote{rational number} (vector), of class \code{"bigq"}.}
  \item{b}{base: from 2 to 36}
  \item{...}{additional arguments passed to methods}
  \item{quote}{(for printing:) logical indicating if the numbers
    should be quoted (as characters are);  the default used to be
    \code{TRUE} (implicitly) till 2011.}
  \item{initLine}{(for printing:) logical indicating if an
    \bold{init}ial line (with the class and length or dimension) should
    be printed. % not yet smart default .. The default prints it for
    % those cases where the class is not easily discernable from the
    % print output.
  }

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

# and, since May 2012,
x ^ 20
stopifnot(is.bigq(x), is.bigq(x + y),
	  x ^ 20 == as.bigz(7)^20 / 2^20)

# convert to string, double
as.character(x)
as.double(x)

# Depict the "S4-class" bigq, i.e., the formal (S4) methods:
if(require("Rmpfr")) # mostly interesting there
  showMethods(class="bigq")
}

\keyword{arith}
