\name{gmp.add}
\alias{gmp.add}
\alias{gmp.sub}
\alias{gmp.mul}
\alias{gmp.divq}
\alias{gmp.divr}
\alias{gmp.div}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Basic arithmetic operators for large integers}
\description{
  Addition, substraction, multiplication, quotient of division,
  remainder of division and division.
  
}
\usage{
gmp.add(a, b)
gmp.sub(a, b)
gmp.mul(a, b)
gmp.divq(a, b)
gmp.divr(a, b)
gmp.div(a, b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{Integer or string from an integer}
  \item{b}{Integer or string from an integer}
}
\details{

}
\value{
  A string of the integer value.
  gmp.div return a float (double).
}

\references{Gnu MP Library see http://swox.com/gmp}

%\references{ ~put references to the literature/web site here ~ }
\author{Antoine Lucas}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
gmp.add(123131312,123131344)
gmp.sub(123131312,123131344)
gmp.mul(123131312,123131344)
gmp.divq(1234,400)
gmp.divr(1234,400)

}
\keyword{arith}% at least one, from doc/KEYWORDS

