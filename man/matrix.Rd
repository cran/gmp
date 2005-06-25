\name{matrix}
\alias{matrix}
\alias{matrix.default}
\alias{matrix.bigz}
\alias{matrix.bigq}
\alias{as.matrix.bigz}
\alias{as.vector.bigq}
\alias{as.vector.bigz}
\alias{\%*\%}
\alias{\%*\%.default}
\alias{\%*\%.bigq}
\alias{\%*\%.bigz}
\alias{ncol.bigq}
\alias{ncol.bigz}
\alias{nrow.bigq}
\alias{nrow.bigz}
\alias{cbind.bigz}
\alias{cbind.bigq}
\alias{rbind.bigz}
\alias{rbind.bigq}
\alias{t.bigq}
\alias{t.bigz}
\alias{dim.bigq}
\alias{dim.bigz}
  
\title{Matrix manipulation with gmp}
\description{
  Overload of all standard tools usefull for matrix manipulation adapted
  to large numbers.
}
\usage{
matrix.bigz(data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL, mod = NA,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{An optional data vector}
  \item{nrow}{the desired number of rows }
  \item{ncol}{the desired number of columns}
  \item{byrow}{logical. If 'FALSE' (the default) the matrix is filled by
          columns, otherwise the matrix is filled by rows.}
  \item{dimnames}{Not implemented}
  \item{mod}{Optional modulus}
  \item{...}{Not used}
}
\details{
  Extract function is the same use for vector or matrix. Then,
  \code{x[i]} return same value as \code{x[i,]}.
  
  Special features concerning bigz class: modulus can be
  \item{Unset}{Just play with large numbers}
  \item{Set with a vector of size 1}{Example:
    \code{matrix.bigz(1:6,nrow=2,ncol=3,mod=7)}
    This mean you work
    in Z/nZ, for the whole matrix. It is the only case
    where \%*\% and solve functions will work in Z/nZ}
  \item{Set with a vector smaller than data}{Example:
    \code{matrix.bigz(1:6,nrow=2,ncol=3,mod=1:5)} then modulus
    is repeated to the end of data. This can be used to define matrix
    with a different modulus at each row
  }
  \item{Set with same size as data}{Modulus is defined for each cell}
}
\value{
  A matrix of class bigz or bigq
}




\author{Antoine Lucas}

\seealso{Solving linear algebra system \code{\link{solve.bigz}} }
\examples{


}

\keyword{arith}

