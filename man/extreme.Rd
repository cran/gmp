\name{Extremes}
\alias{max.bigz}
\alias{max.bigq}
\alias{min.bigz}
\alias{min.bigq}


\title{Returns the maxima and minima of the input values.}


\usage{
     max.bigz(..., na.rm=FALSE)
     max.bigq(..., na.rm=FALSE)
     min.bigz(..., na.rm=FALSE)
     min.bigq(..., na.rm=FALSE)
}

\description{
  Theses functions overload min and max function for big rationals and
  big integers
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
  x = as.bigz(1:10)
  ##_call max.bigz => return 10
  max(x)
  min(x)

  x = as.bigq(1:10)
  ##_call max.bigq
  max(x)
  min(x)
  
}

\seealso{\code{\link[base]{max}}}

\keyword{arith}

