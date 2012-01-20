\name{Extremes}
\alias{max.bigz}
\alias{max.bigq}
\alias{min.bigz}
\alias{min.bigq}

\title{Extrema (Maxima and Minima)}

\usage{
\method{max}{bigz}(..., na.rm=FALSE)
\method{max}{bigq}(..., na.rm=FALSE)
\method{min}{bigz}(..., na.rm=FALSE)
\method{min}{bigq}(..., na.rm=FALSE)
}

\description{
  We provide S3 \code{\link{methods}} for \code{\link{min}} and
  \code{\link{max}} for big rationals (\code{bigq}) and big integers
  (\code{biqz});  consequently, \code{\link{range}()} works as well.
}

\arguments{
  \item{...}{numeric arguments}
  \item{na.rm}{a logical indicating whether missing values should be removed.}
}

\value{
  return an element of class bigz or bigq.
}

\author{Antoine Lucas}

\examples{
 x <- as.bigz(1:10)
 max(x)
 min(x)
 range(x) # works correctly via default method

 Q <- as.bigq(1:10, 3)
 max(Q)
 min(Q)

stopifnot(range(x) == c(1,10), 3*range(Q) == c(1,10))
}

\seealso{\code{\link[base]{max}}}

\keyword{arith}

