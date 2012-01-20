\name{bigz operators}
\alias{add.bigz}
\alias{+.bigz}
\alias{sub.bigz}
\alias{-.bigz}
\alias{mul.bigz}
\alias{*.bigz}
\alias{div.bigz}
\alias{/.bigz}
\alias{mod.bigz}
\alias{\%\%.bigz}
\alias{\%/\%.bigz}
\alias{divq.bigz}
\alias{abs.bigz}
\alias{inv.bigz}
\alias{inv}
\alias{pow.bigz}
\alias{pow}
\alias{^.bigz}
\alias{log.bigz}
\alias{log2.bigz}
\alias{log10.bigz}
\title{Basic Arithmetic Operators for Large Integers ("bigz")}
\description{
  Addition, substraction, multiplication, (integer) division,
  remainder of division, multiplicative inverse,
  power and logarithm functions.
}
\usage{
add.bigz(e1, e2)
sub.bigz(e1, e2 = NULL)
mul.bigz(e1, e2)
div.bigz(e1, e2)
divq.bigz(e1,e2) ## ==  e1 \%/\% e2
mod.bigz(e1, e2) ## ==  e1 \%\%  e2
\method{abs}{bigz}(x)

inv.bigz(a, b,...)## == (1 / a) (modulo b)
% inv(a, ...)
pow.bigz(e1, e2,...)## == e1 ^ e2
% pow(a, ...)
\S3method{log}{bigz}(x, base=exp(1))
\S3method{log2}{bigz}(x)
\S3method{log10}{bigz}(x)
}

\arguments{
  \item{x}{bigz, integer or string from an integer}
  \item{e1, e2, a,b}{bigz, integer or string from an integer}
  \item{base}{base of the logarithm; base e as default}
  \item{...}{Additional parameters}
}
\details{
  For details about the internal modulus state, see the manpage of
  \code{\link{bigz}}.

  div or "/" return a rational number; divq or "\%/\%" return the
  quotient of integer division.

  Operators can be used directly when objects are of class bigz:
  a + b, log(a), etc.
}
\value{
  A bigz class representing the result of the arithmetic operation.
}

\references{The GNU MP Library, see \url{http://gmplib.org}}

\author{Immanuel Scholz and Antoine Lucas}

\examples{
# 1+1=2
as.bigz(1) + 1
as.bigz(2)^10
as.bigz(2)^200

# if my.large.num.string is set to a number, this returns the least byte
(my.large.num.string <- paste(sample(0:9, 200, replace=TRUE), collapse=""))
mod.bigz(as.bigz(my.large.num.string), "0xff")

# power exponents can be up to MAX_INT in size, or unlimited if a
# bigz's modulus is set.
pow.bigz(10,10000)

## Modulo 11,   7 and 8 are inverses :
as.bigz(7, mod = 11) * 8 ## ==>  1  (mod 11)
inv.bigz(7, 11)## hence, 8
a <- 1:10
(i.a <- inv.bigz(a, 11))
d <- as.bigz(7)
a \%/\% d  # = divq(a, d)
a \%\%  d  # = mod.bigz (a, d)

stopifnot(inv.bigz(7, 11) == 8,
          all(as.bigz(i.a, 11) * a == 1),
          identical(a \%/\% d, divq.bigz(1:10, 7)),
          identical(a \%\%  d, mod.bigz (a, d))
 )

}

\keyword{arith}
%\keyword{bigz}
