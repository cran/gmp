\name{add.bigz}
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
a + b (if a or b is a bigz vector)
sub.bigz(a, b)
a - b (if a or b is a bigz vector)
mul.bigz(a, b)
a * b (if a or b is a bigz vector)
div.bigz(a, b)
a / b (if a or b is a bigz vector)
mod.bigz(a, b)
a \%\% b (if a or b is a bigz vector)
inv.bigz(a, b,...)
inv(a, ...)
pow.bigz(a, b,...)
pow(a, ...)
a ^ b (if a or b is a bigz vector)
log.bigz(x, base=exp(1))
log(a, base=exp(1)) (if a is a bigz vector)
log2.bigz(a)
log2(a) (if a is a bigz vector)
log10.bigz(a)
log10(a) (if a is a bigz vector)
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
}
\value{
  A bigz class representing the result of the arithmetic operation.
}


\references{
   \item Gnu MP Library see \url{http://swox.com/gmp},
   \item Home page: \url{http://mulcyber.toulouse.inra.fr/projects/gmp/}
 }


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

