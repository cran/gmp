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


as.matrix.bigq <- function(x, ...)
  {
    n <- length(x)
    p <- 1
    if((class(x) == ("matrix") )| (class(x) == "data.frame") | (class(x) == "bigq"))
      {
        n=dim(x)[1]
        p=dim(x)[2]        
      }
    matrix.bigz(x,nrow=n,ncol=p)
  }


as.vector.bigq <- function(x, mode="any")
  {
    attr(x,"nrow")<- NULL
    return(x)
  }


dim.bigq <- function(x)
  return(c(attr(x,"nrow"),length(x)/attr(x,"nrow")))

"dim<-.bigq" <- function(x,value)
{
  ## TODO: check...
  attr(x,"nrow") <- as.integer(value[1])
  x
}

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



apply.bigq <- function(X, MARGIN,FUN)
{
  ## This change matrix to a list.
  X = .Call("gmpMatToListQ",X,as.integer(MARGIN),PACKAGE="gmp")
  ## here: std lapply
  lst = lapply(X,FUN)

  ## to change list to vector
  .Call("bigrational_c",lst)
}


