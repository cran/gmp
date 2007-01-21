\name{bigz}
\alias{bigz}
\alias{as.bigz}
\alias{as.character.bigz}
\alias{as.double.bigz}
\alias{print.bigz}
\alias{is.na.bigz}
\title{Large sized integer values}
\description{
  Type class supporting arithmetic operations on very large integers.
}
\usage{
as.bigz(a, mod = NA)
as.character.bigz(x,b=10,...)
as.double.bigz(x,...)
%as.bigz(str, mod = NA)
%as.bigz(big, mod = NA)
is.na.bigz(x)
print.bigz(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%  \item{int}{Integer value (or logical, which will be interpreted as 1 resp. 0)}
%  \item{num}{Numeric value}
%  \item{str}{String value, either starting with \code{0x} for hexadecimal, \code{0b} for
%  binary or without prefix for decimal values. Any format error results in \code{0}}
%  \item{big}{A bigz value (class of type \code{"\link{bigz}"})}
\item{a}{Either integer, numeric or string value
  (String value: ither starting with \code{0x} for hexadecimal, \code{0b} for
  binary, \code{0} for octal or without prefix for decimal values.
  Any format error results return an Error message)}
\item{b}{Base: from 2 to 32}
\item{x}{Numeric value}
\item{...}{Additional parameters}
\item{mod}{An integer, numeric, string or bigz of the internal modulus, see below.}
}
\value{
  A bigz class representing the parameter value.
}
\details{
  Bigzs are integers of infinite, but given length (means: only
  restricted by the host memory). Basic arithmetic operations can be performed
  on bigzs as addition, subtraction, multiplication, division, modulation
  (remainder of division), power, multiplicative inverse, calculating of the
  greatest common divisor, test whether the integer is prime and other things
  that comes in need when performing standard cryptographic operations.

  For a review of basic arithmetics, see \code{"\link{add.bigz}"}.

  The most important logical operators are supported, such as \code{"=="}, \code{"!="},
  \code{"<"}, \code{"<="}, \code{">"}, and \code{">="}.

  Objects of class \code{"\link{bigz}"} may have an attribute \code{mod} 
  which specifies a modulus that is applied after each arithmetic operation. 
  If the result is going to have a modulus, 
  
    \code{result = mod.bigz(result, modulus)}
    
  is called after performing the arithmetic operation and the result will have the 
  attribute \code{mod} set accordingly.

  Powers of bigzs can only be performed, if either a modulus is going to
  be applied to the result bigz or if the exponent fits into an integer
  value. So if you want to calculate a power in a finite group, don't use 
  \code{a ^ b \%\% c}, but use \code{as.bigz(a,c) ^ b}, instead.
  
  The following rules for the result's modulus apply when performing
  arithmetic operations on bigzs:

  \item If none of the operand has a modulus set, the result will not have a modulus.
  \item If both operands have a different modulus, the result will not have a
  modulus, except in case of \code{"\link{mod.bigz}"}, where the second operands 
  value is used.
  \item If only one of the operands has a modulus or both have a common (the
  same), it is set and applied to the result, except in case of mod.bigz,
  where the second operands value is used.
}

\note{
  x = as.bigz(1234567890123456789012345678901234567890)
  will not work as integer will be first converted to a double
  and then converted to a "bigz" element.

  Correct syntaxe:
  x = as.bigz("1234567890123456789012345678901234567890")
}


\references{Gnu MP Library see http://swox.com/gmp}

%\references{ ~put references to the literature/web site here ~ }
\author{Immanuel Scholz}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## 1+1=2
a = as.bigz(1)
a + a

\dontrun{
## calculate c = x^e mod n
x <- as.bigz("0x123456789abcdef") # my secret message
e <- as.bigz(3) # something smelling like a dangerous public RSA exponent
n <- as.bigz("0x4378a27...")  # probably a product of two primes

# first way to do it right
modulus(x) <- n
c <- x ^ e

# similar second way (maybe more sensefull if you reuse e) to do it right
modulus(e) <- n
c <- x ^ e

# third way to do it right
c <- x ^ as.bigz(e, n)

# WRONG! (although very beautiful. Maybe ok for small examples)
c <- x ^ e %% n

# Return result in hexa
as.character(c,b=16)
}
}

\keyword{arith}
%\keyword{bigz}

