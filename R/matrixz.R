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
    if((class(x) == ("matrix") )| (class(x) == "data.frame") | (class(x) == "bigz"))
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

  
##aperm.bigz <- function(a,perm, resize= TRUE)
##  {
##    dims <- dim(a)
##    if (missing(perm)) 
##      perm <- c(1,2)
##    if(perm[1] > perm[2])
##      ans = .Call("bigint_transposeR",
##        a,    
##       PACKAGE="gmp")
##   else
##     ans = a  
##    if(!resize)
##      dim(ans) <- dims
##    ans
##}


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

"dim<-.bigz" <- function(x,value)
{
  ## TODO: check...
  attr(x,"nrow") <- as.integer(value[1])
  x
}

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

apply <- function(X, MARGIN,FUN)
  UseMethod("apply")

apply.default <- function(X, MARGIN,FUN)
  base::apply(X, MARGIN,FUN)


apply.bigz <- function(X, MARGIN,FUN)
{
  ## This change matrix to a list.
  X = .Call("gmpMatToListZ",X,as.integer(MARGIN),PACKAGE="gmp")
  ## here: std lapply
  lst = lapply(X,FUN)

  ## to change list to vector
  .Call("biginteger_c",lst)
}


