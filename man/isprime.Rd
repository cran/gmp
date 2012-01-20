\name{isprime}
\alias{isprime}
\title{Determine if number is (very probably) prime}
\description{
  Determine whether the number is prime or not. 2: number is prime,
  1, number is probably prime (without beeing certain), 0 number is
  composite.
}
\usage{
isprime(n, reps = 40)  % we feel better with 40 :)
}
\arguments{
  \item{n}{integer number, to be tested}
  \item{reps}{integer number of repeats}
}
\details{
  This function does some trial divisions, then some Miller-Rabin
  probabilistic primary tests.  \code{reps} controls how many such tests are
  done, 5 to 10 is a resonable number.  More will reduce the chances
  of a composite being returned as \dQuote{probably prime}.

}
\value{
  \item{0}{\eqn{n} is not prime}
  \item{1}{\eqn{n} is probably prime}
  \item{2}{\eqn{n} is prime}
}
\references{The GNU MP Library, see \url{http://gmplib.org}}

\author{Antoine Lucas}

\seealso{\code{\link{nextprime}}
}
\examples{
isprime(210)
isprime(71)

# All primes numbers from 1 to 100
t <- isprime(1:99)
(1:99)[t > 0]

table(isprime(1:10000))# 0 and 2 : surely prime or not prime

primes <- function(n) {
  ## all primes <= n
  stopifnot(length(n) == 1, n <= 1e7) # be reasonable
  p <- c(2L, as.integer(seq(3, n, by=2)))
  p[isprime(p) > 0]
}

## quite quickly, but for these small numbers
## still slower than e.g., sfsmisc::primes()
system.time(p100k <- primes(100000))

## The first couple of Mersenne primes:
p.exp <- primes(1000)
Mers <- as.bigz(2) ^ p.exp - 1
isp.M <- sapply(seq_along(Mers), function(i) isprime(Mers[i], reps=256))
cbind(p.exp, isp.M)[isp.M > 0,]
Mers[isp.M > 0]
}
\keyword{arith}

