\name{Bigq operators}
\alias{add.bigq}
\alias{+.bigq}
\alias{sub.bigq}
\alias{-.bigq}
\alias{mul.bigq}
\alias{*.bigq}
\alias{div.bigq}
\alias{/.bigq}
\alias{\%/\%.bigq}
\alias{abs.bigq}

\title{Basic arithmetic operators for large rationals}
\description{
  Addition, substraction, multiplication, division.

}
\usage{
add.bigq(a, b)
sub.bigq(a, b=NULL)
mul.bigq(a, b)
div.bigq(a, b)
abs.bigq(a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a,b}{bigq, integer or string from an integer}
  \item{...}{Additional parameters}
}
\details{
  Operators can be use directly when objects are of class bigq:
  a + b, a * b, etc.

}
\value{
  A bigq class representing the result of the arithmetic operation.
}

\references{http://mulcyber.toulouse.inra.fr/projects/gmp/}

%\references{ ~put references to the literature/web site here ~ }
\author{Immanuel Scholz and Antoine Lucas}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
# 1/3+1=4/3
as.bigq(1,3) + 1


}

\keyword{arith}
%\keyword{bigz}

