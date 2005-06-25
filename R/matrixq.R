matrix.bigq <- function(data=NA,nrow=1, ncol=1, byrow=FALSE,dimnames =NULL,
                        den=NA,...)
  {
    .Call("as_matrixq",data,
          as.integer(nrow),
          as.integer(ncol),
          as.integer(byrow),
          den,
          PACKAGE="gmp")
  }


as.matrix.bigz <- function(x)
  {
    .Call("as_matrixq",x,
          as.integer(1),
          as.integer(1),
          FALSE,
          NA,
          PACKAGE="gmp")
  }


as.vector.bigq <- function(x, mode="any")
  {
    attr(x,"nrow")<- NULL
    return(x)
  }


dim.bigq <- function(x)
  return(c(attr(x,"nrow"),length(x)/attr(x,"nrow")))

nrow.bigq <- function(x)
  return(attr(x,"nrow"))

ncol.bigq <- function(x)
  return(length(x)/attr(x,"nrow"))

t.bigq <- function(x)
  {
    .Call("bigq_transposeR",
          x,    
          PACKAGE="gmp")
  }


"%*%.bigq" <- function(x,y)
  {
    .Call("matrix_mul_q",
          x,
          y,
          PACKAGE="gmp")
  }


cbind.bigq <- function(..., recursive = FALSE)
  {
    a <- list(...)

    return(.Call("bigrational_cbind",a,PACKAGE = "gmp"))
  }

rbind.bigq <- function(..., recursive = FALSE)
  {
    a <- list(...)

    return(.Call("bigrational_rbind",a,PACKAGE = "gmp"))
  }
