#----------------------------------------------------------
#
#  Author        : Immanuel Scholz (immanuel.scholz@gmx.de)
#		   Technische Universitaet Dresden
#
#  Brief         : Stub to call the dll functions
#
#  Licence       : GPL 
#
#----------------------------------------------------------


"+.bigz" <- function(...) add.bigz(...)
add.bigz<- function(a,b)
{
    .Call("biginteger_add",a,b, PACKAGE= "gmp")
}

"-.bigz" <- function(...) sub.bigz(...)
sub.bigz<- function(a,b=NULL)
{
    if(is.null(b))
      .Call("biginteger_sub",0,a, PACKAGE= "gmp")
    else
      .Call("biginteger_sub",a,b, PACKAGE= "gmp")
}

"*.bigz" <- function(...) mul.bigz(...)
mul.bigz<- function(a,b)
{
    .Call("biginteger_mul",a,b, PACKAGE= "gmp")
}

"%/%.bigz" <- function(...) divq.bigz(...)
divq.bigz<- function(a,b)
{
    .Call("biginteger_div",a,b, PACKAGE= "gmp")
}

"/.bigz" <- function(...) div.bigz(...)
div.bigz<- function(a,b)
{
  ismod <- FALSE
  if(class(a) == "bigz")
    if(!is.null(modulus(a)) )
      ismod =TRUE
  if(class(b) == "bigz")
    if(!is.null(modulus(b)) )
      ismod =TRUE
  if(ismod)
    .Call("biginteger_div",a,b, PACKAGE= "gmp")
  else
    .Call("bigrational_as",a,b, PACKAGE= "gmp")
}

"%%.bigz" <- function(...) mod.bigz(...)
mod.bigz<- function(a,b)
{
    .Call("biginteger_mod",a,b, PACKAGE= "gmp")
}

"^.bigz" <- function(...) pow.bigz(...)
pow.bigz<- function(a,b,...)
{
    .Call("biginteger_pow",a,b, PACKAGE= "gmp")
}

inv.bigz<- function(a,b,...)
{
    .Call("biginteger_inv",a,b, PACKAGE= "gmp")
}


gcd <- function(a,b)
      UseMethod("gcd")

gcd.default <- function(a,b)
  return(as.integer(gcd.bigz(a,b)))

gcd.bigz<- function(a,b)
{
    .Call("biginteger_gcd",a,b, PACKAGE="gmp")
}

##lcm <- function(a,b)
##      UseMethod("lcm")

lcm.default <- function(a,b)
  return(as.integer(lcm.bigz(a,b)))

lcm.bigz<- function(a,b)
{
    .Call("biginteger_lcm",a,b, PACKAGE="gmp")
}

print.bigz<- function(x,...)
{
  if (length(x)>0)
    print(as.character(x))
  else
    cat("bigz(0)\n")
}

as.bigz<- function(a,mod = NA)
{
  if(class(a) == "bigq")
    as.bigz.bigq(a,mod)
  else
  .Call("biginteger_as", a, mod, PACKAGE="gmp")
}

as.character.bigz<- function(x,b=10,...)
{
    .Call("biginteger_as_character", x,as.integer(b), PACKAGE="gmp")
}

as.double.bigz<- function(x,...)
{
    .Call("biginteger_as_numeric", x, PACKAGE="gmp")
}

"[[.bigz"<- function(a,b=NA)
{
    .Call("biginteger_get_at", a, b, PACKAGE="gmp")
}

"[[<-.bigz"<- function(dst, idx=NA, value)
{
    .Call("biginteger_set_at", dst, idx, value, PACKAGE="gmp")
}

length.bigz<- function(x)
{
    .Call("biginteger_length", x, PACKAGE="gmp")
}

"length<-.bigz"<- function(x, value)
{
    .Call("biginteger_setlength", x, value, PACKAGE="gmp")
}

modulus <- function(a, value) 
{
    UseMethod("modulus")
}

"modulus<-" <- function (a, value) 
{
    UseMethod("modulus<-")
}

inv <- function(a,...) 
{
    UseMethod("inv")
}

pow <- function(a,...) 
{ 
    UseMethod("pow")
}

modulus.bigz <- function(a,value) 
{
    attr(a, "mod")
}

"modulus<-.bigz" <- function(a, value) 
{
    as.bigz(a,value)
}

is.na.bigz <- function(x) 
{
    .Call("biginteger_is_na", x, PACKAGE="gmp")
}

"<.bigz" <- function(a,b) 
{
    .Call("biginteger_lt", a, b, PACKAGE="gmp")
}

">.bigz" <- function(a,b) 
{
    .Call("biginteger_gt", a, b, PACKAGE="gmp")
}

"<=.bigz" <- function(a,b) 
{
    .Call("biginteger_lte", a, b, PACKAGE="gmp")
}

">=.bigz" <- function(a,b) 
{
    .Call("biginteger_gte", a, b, PACKAGE="gmp")
}

"==.bigz" <- function(a,b) 
{
    .Call("biginteger_eq", a, b, PACKAGE="gmp")
}

"!=.bigz" <- function(a,b) 
{
    .Call("biginteger_neq", a, b, PACKAGE="gmp")
}

abs.bigz <- function(a)
  {
    .Call("biginteger_abs",a,PACKAGE="gmp")
  }

sign.bigz <- function(a)
  {
    .Call("biginteger_sgn",a,PACKAGE="gmp")
  }


max.bigz <- function(...,na.rm=FALSE)
{
 .Call("biginteger_max",c.bigz(...),na.rm,PACKAGE="gmp")
}

min.bigz <- function(...,na.rm=FALSE)
{
 .Call("biginteger_min",c.bigz(...),na.rm,PACKAGE="gmp")
}

cumsum.bigz <- function(...)
{
 .Call("biginteger_cumsum",c.bigz(...),PACKAGE="gmp")
}

prod <- function(...)
{
  UseMethod("prod")
}

prod.default <- function(...)
{
  base::prod(...)
}

prod.bigz <- function(...)
{
 .Call("biginteger_prod",c.bigz(...),PACKAGE="gmp")
}



# TODO: Compute log functions directly on bigzs to increase accuracy

log.bigz <- function(x,base=exp(1))
{
     log(as.double(x), base)
}

log2.bigz <- function(a)
{
     log(as.double(a), 2)
}

log10.bigz <- function(a)
{
     log(as.double(a), 10)
}

c.bigz <- function(..., recursive = FALSE)
  {
    a <- list(...)
    .Call("biginteger_c",a,PACKAGE = "gmp")
  }
rep.bigz <- function(x,times,...)
  {

    .Call("biginteger_rep",x,times,PACKAGE = "gmp")
  }


# Isprime, return:
#   0 if not prime
#   1 if probably prime
#   2 if prime
isprime <- function(n,reps=40)
  {
    .Call("biginteger_is_prime",
          n,
          as.integer(reps),
          PACKAGE="gmp")
  }

nextprime<- function(n)
 {
   .Call("biginteger_nextprime",n,
         PACKAGE= "gmp"
      ) 
 }

gcdex <- function(a,b)
 {
   .Call("biginteger_gcdex",a,b,
         PACKAGE= "gmp"
      ) 
 }

urand.bigz<- function(nb=1,size=200,seed=0)
  {
    if(seed==0)
      ok <- 0
    else
      ok <- 1

    .Call("biginteger_rand_u",
          as.integer(nb),
          as.integer(size),
          seed,
          as.integer(ok),
       PACKAGE= "gmp"
       )
  }

sizeinbase <- function(a,b=10)
 {
   if(as.integer(b)<2)
     {
       print("base must be >= 2")
       return(0)
     }
   .Call("biginteger_sizeinbase",a,as.integer(b),
         PACKAGE= "gmp"
      ) 
 }


fibnum <- function (n)
{
   .Call("fibnum",as.integer(n),
         PACKAGE= "gmp"
      ) 

}
fibnum2 <- function (n)
{
   .Call("fibnum2",as.integer(n),
         PACKAGE= "gmp"
      ) 

}
 lucnum <- function (n)
{
   .Call("lucnum",as.integer(n),
         PACKAGE= "gmp"
      ) 

}
lucnum2 <- function (n)
{
   .Call("lucnum2",as.integer(n),
         PACKAGE= "gmp"
      ) 

}
  
factorize <- function(n) 
{
   .Call("factorR",as.bigz(n),
         PACKAGE= "gmp"
      ) 

}

## overload as.vector
as.vector.bigz <- function(x,mode="any")
  return(x)



solve.bigz <- function(a,b,...)
  {
    if(missing(b))
      .Call("inverse_z",a,PACKAGE="gmp")
    else
      .Call("solve_z",a,b,PACKAGE="gmp")
  }


"[.bigz"<- function(a,b=NULL,c=NULL)
{
  .Call("matrix_get_at_z", a, b,c, PACKAGE="gmp")
}

"[<-.bigz"<- function(dst,idx=NULL,jdx=NULL,value)
{
  .Call("matrix_set_at_z", dst, value,idx,jdx , PACKAGE="gmp")
}
