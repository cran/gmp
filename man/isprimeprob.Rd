\name{isprimeprob}
\alias{isprimeprob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determine if number is prime}
\description{
  Determine whether the number is prime or not. 2: number is prime,
  1, number is probably prime (without beeing certain), 0 number is
  composite.
}
\usage{
isprimeprob(n, reps = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Integer number, to be tested}
  \item{reps}{Integer, number of repeats }
}
\details{
  This function does some trial divisions, then some Miller-Rabin
  probabilisticprimary tests. reps controls how many such tests are
  done, 5 to 10 is a resonable number. More will reduce the chances
  of a composite being returned as "probably prime". 

}
\value{
  \item{0}{Number is not prime}
  \item{1}{Number is probably prime}
  \item{2}{Number is prime}


}
\references{Gnu MP Library see http://swox.com/gmp}

\author{Antoine Lucas}



\seealso{\code{\link{nextprime}}}
\examples{
isprimeprob(210)
isprimeprob(71)

# All primes numbers from 1 to 100
t <-sapply(1:100,isprimeprob)
(1:100)[t>0]


}
\keyword{arith }% at least one, from doc/KEYWORDS

