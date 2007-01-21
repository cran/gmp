\name{apply}
\alias{apply}
\alias{apply.bigz}
\alias{apply.bigq}
\alias{apply.default}
\title{apply over a matrix}

\usage{
apply(X, MARGIN, FUN)
}


\description{
  Theses functions overload apply function for big rationals and
  big integers
}

\arguments{
  \item{X}{A matrix of class bigz or bigq}
  \item{MARGIN}{1: apply function on rows; 2: apply function of columns}
  \item{FUN}{function to be applied}
}


\value{
  return a vector of class bigz or bigq.
}

\author{Antoine Lucas}

\examples{
  x = as.bigz(matrix(1:12,3))
  apply(x,1,min)
  apply(x,2,max)
  
  x = as.bigq(matrix(1:12,3))
  apply(x,1,min)
  apply(x,2,max)
  
}

\keyword{arith}


\seealso{\code{\link[base]{apply}}}

