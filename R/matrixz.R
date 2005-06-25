matrix <- function(...)
  UseMethod("matrix")

matrix.default <- function(...)
{
  base::matrix(...)
}

matrix.bigz <- function(data=NA,nrow=1, ncol=1, byrow=FALSE,dimnames =NULL,
                        mod=NA,...)
  {
    .Call("as_matrixz",data,
          as.integer(nrow),
          as.integer(ncol),
          as.integer(byrow),
          mod,
          PACKAGE="gmp")
  }


as.matrix.bigz <- function(x)
  {
    n <- length(x)
    p <- 1
    if((class(x) == ("matrix") )| (class(x) == "data.frame") | (class(x) == "matrixz"))
      {
        n=dim(x)[1]
        p=dim(x)[2]        
      }
    matrix.bigz(x,nrow=n,ncol=p)
  }


as.vector.bigz <- function(x, mode="any")
  {
    attr(x,"nrow")<- NULL
    return(x)
  }

t.bigz <- function(x)
  {
    .Call("bigint_transposeR",
          x,    
          PACKAGE="gmp")
  }

"%*%" <- function(x,y)
  UseMethod("%*%")

"%*%.default" <- function(x,y)
{
  base::"%*%"(x,y)
}

"%*%.bigz" <- function(x,y)
  {
    .Call("matrix_mul_z",
          x,
          y,
          PACKAGE="gmp")
  }

dim.bigz <- function(x)
  return(c(attr(x,"nrow"),length(x)/attr(x,"nrow")))

nrow.bigz <- function(x)
  return(attr(x,"nrow"))

ncol.bigz <- function(x)
  return(length(x)/attr(x,"nrow"))


cbind.bigz <- function(..., recursive = FALSE)
  {
    a <- list(...)
    return(.Call("biginteger_cbind",a,PACKAGE = "gmp"))
  }

rbind.bigz <- function(..., recursive = FALSE)
  {
    a <- list(...)
    return(.Call("biginteger_rbind",a,PACKAGE = "gmp"))
  }
