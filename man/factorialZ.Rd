\name{factorialZ}
\alias{factorialZ}
\alias{chooseZ}
\title{Factorial and Binomial Coefficient as Big Integer}
\description{
  Efficiently compute the factorial \eqn{n!} or a binomial coefficient
  \eqn{{n\choose k}}{choose(n, k)} as big integer (class \code{\link{bigz}}).
}
\usage{
factorialZ(n)
chooseZ(n, k)
}
\arguments{
  \item{n}{non-negative integer (scalar for now), for
    \code{factorialZ}.  For \code{chooseZ}, may be a \code{bigz} big
    integer, also negative.}
  \item{k}{non-negative integer.}
}
\value{
  a number of (S3) class \code{\link{bigz}}.
}
\seealso{
  \code{\link{factorial}} and \code{\link{gamma}} in base \R;
}
\examples{
factorialZ(200)

n <- 1000
f1000 <- factorialZ(n)
stopifnot(1e-15 > abs(as.numeric(1 - lfactorial(n)/log(f1000))))

system.time(replicate(8, f1e4 <<- factorialZ(10000)))
nchar(as.character(f1e4))# 35660 ... (too many to even look at ..)

chooseZ(1000, 100)
chooseZ(as.bigz(2)^120, 10)
}
\keyword{arith}
