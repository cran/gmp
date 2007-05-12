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
\alias{divq.bigz}
\alias{\%/\%.bigz}
\alias{abs.bigz}
\alias{inv.bigz}
\alias{inv}
\alias{pow.bigz}
\alias{pow}
\alias{^.bigz}
\alias{log.bigz}
\alias{log2.bigz}
\alias{log10.bigz}
\title{Basic arithmetic operators for large integers}
\description{
  Addition, substraction, multiplication, division,
  remainder of division, multiplicative inverse, 
  power and logarithm functions.
}
\usage{
add.bigz(a, b)
sub.bigz(a, b = NULL)
mul.bigz(a, b)
div.bigz(a, b)
divq.bigz(a,b)
mod.bigz(a, b)
abs.bigz(x)
inv.bigz(a, b,...)
inv(a, ...)
pow.bigz(a, b,...)
pow(a, ...)
log.bigz(x, base=exp(1))
log2.bigz(a)
log10.bigz(a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{bigz, integer or string from an integer}
  \item{a}{bigz, integer or string from an integer}
  \item{b}{bigz, integer or string from an integer}
  \item{base}{base of the logarithm; base e as default}
  \item{...}{Additional parameters}
}
\details{
  For details about the internal modulus state, see the manpage of
  "bigz".

  div or "/" return a rational number; divq or "\%/\%" return the
  quotient of division.

  Operators can be use directly when objects are of class bigz:
  a + b, log(a), etc.
}
\value{
  A bigz class representing the result of the arithmetic operation.
}

\references{Gnu MP Library see http://swox.com/gmp}

%\references{ ~put references to the literature/web site here ~ }
\author{Immanuel Scholz and Antoine Lucas}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
# 1+1=2
as.bigz(1) + 1

# if my_large_number_string is set to a number, it returns the least byte
\dontrun{
mod.bigz(as.bigz(my_large_number_string),"0xff")
}

# power exponents can be up to MAX_INT in size, or unlimited if a 
# bigz's modulus is set.
pow.bigz(10,10000)
}

\keyword{arith}
%\keyword{bigz}

