#----------------------------------------------------------
#
#  Author        : Antoine Lucas (adapted from biginteger class made by
#                                 Immanuel Scholz)
#
#  Brief         : Stub to call the dll functions
#
#  Licence       : GPL 
#
#----------------------------------------------------------


"+.bigq" <- function(...) add.bigq(...)
add.bigq<- function(a,b)
{
    .Call("bigrational_add",a,b, PACKAGE= "gmp")
}

"-.bigq" <- function(...) sub.bigq(...)
sub.bigq<- function(a,b=NULL)
{
  if(is.null(b))
    .Call("bigrational_sub",0,a, PACKAGE= "gmp")
  else
    .Call("bigrational_sub",a,b, PACKAGE= "gmp")
}

"*.bigq" <- function(...) mul.bigq(...)
mul.bigq<- function(a,b)
{
    .Call("bigrational_mul",a,b, PACKAGE= "gmp")
}

"/.bigq" <- function(...) div.bigq(...)
div.bigq<- function(a,b)
{
    .Call("bigrational_div",a,b, PACKAGE= "gmp")
}


print.bigq<- function(x,...)
{
  if (length(x)>0)
    print(as.character(x))
  else
    cat("bigrational(0)\n")
}

as.bigq<- function(n,d=1)
{
    .Call("bigrational_as", n, d, PACKAGE="gmp")
}

as.character.bigq<- function(a,b=10)
{
    .Call("bigrational_as_character", a,as.integer(b), PACKAGE="gmp")
}

as.double.bigq<- function(x,...)
{
    .Call("bigrational_as_numeric", x, PACKAGE="gmp")
}


denominator <- function(x)
  .Call("bigrational_den",x,PACKAGE="gmp")

"denominator<-" <- function(x,value)
  as.bigq(numerator(x),value)


numerator <- function(x)
  .Call("bigrational_num",x,PACKAGE="gmp")

"numerator<-" <- function(x,value)
  as.bigq(value,denominator(x))


as.bigz.bigq<- function(a,mod = NA)
{  
  as.bigz(numerator(a) %/% denominator(a),mod)
}

abs.bigq <- function(a)
   {
     numerator(a) <- abs(numerator(a))
     a
   }

sign.bigq <- function(a)  
  sign(numerator(a))
  


"[.bigq"<- function(a,b=NA)
{
    .Call("bigrational_get_at", a, b, PACKAGE="gmp")
}

"[<-.bigq"<- function(dst, idx=NA, value)
{
    .Call("bigrational_set_at", dst, idx, value, PACKAGE="gmp")
}

length.bigq<- function(a)
{
    .Call("bigrational_length", a, PACKAGE="gmp")
}

"length<-.bigq"<- function(vec, value)
{
    .Call("bigrational_setlength", vec, value, PACKAGE="gmp")
}

is.na.bigq <- function(a) 
{
    .Call("bigrational_is_na", a, PACKAGE="gmp")
}

"<.bigq" <- function(a,b) 
{
    .Call("bigrational_lt", a, b, PACKAGE="gmp")
}

">.bigq" <- function(a,b) 
{
    .Call("bigrational_gt", a, b, PACKAGE="gmp")
}

"<=.bigq" <- function(a,b) 
{
    .Call("bigrational_lte", a, b, PACKAGE="gmp")
}

">=.bigq" <- function(a,b) 
{
    .Call("bigrational_gte", a, b, PACKAGE="gmp")
}

"==.bigq" <- function(a,b) 
{
    .Call("bigrational_eq", a, b, PACKAGE="gmp")
}

"!=.bigq" <- function(a,b) 
{
    .Call("bigrational_neq", a, b, PACKAGE="gmp")
}

c.bigq <- function(..., recursive = FALSE)
  {
    a <- list(...)
    .Call("bigrational_c",a,PACKAGE = "gmp")
  }
rep.bigq <- function(x,times,...)
  {

    .Call("bigrational_rep",x,times,PACKAGE = "gmp")
  }
