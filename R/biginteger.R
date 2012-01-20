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


"+.bigz" <- add.bigz <- function(e1, e2) .Call(biginteger_add, e1, e2)

"-.bigz" <- sub.bigz <- function(e1, e2=NULL)
{
    if(is.null(e2))
      .Call(biginteger_sub, 0, e1)
    else
      .Call(biginteger_sub, e1, e2)
}

"*.bigz" <- mul.bigz <- function(e1, e2) .Call(biginteger_mul, e1, e2)

"%/%.bigz" <- divq.bigz <- function(e1, e2) .Call(biginteger_div, e1, e2)

"/.bigz" <- div.bigz <- function(e1, e2)
{
  ismod <- ((inherits(e1, "bigz") && !is.null(modulus(e1))) ||
            (inherits(e2, "bigz") && !is.null(modulus(e2))))
  if(ismod)
    .Call(biginteger_div, e1, e2)
  else
    .Call(bigrational_as, e1, e2)
}

"%%.bigz" <- mod.bigz <- function(e1, e2) .Call(biginteger_mod,e1, e2)

"^.bigz" <- pow.bigz <- function(e1, e2,...) .Call(biginteger_pow,e1, e2)

inv.bigz <- function(a,b,...) .Call(biginteger_inv,a,b)

gcd <- function(a,b)
      UseMethod("gcd")
gcd.default <- function(a,b) as.integer(gcd.bigz(a,b))
gcd.bigz <- function(a,b) .Call(biginteger_gcd,a,b)

## just because lcm() is a trivial function in 'graphics' .. hmm
##lcm <- function(a,b)
##      UseMethod("lcm")

lcm.default <- function(a,b)
  as.integer(lcm.bigz(a,b))
lcm.bigz <- function(a,b) .Call(biginteger_lcm,a,b)

print.bigz <- function(x, quote = FALSE, ...)
{
  if((n <- length(x)) > 0) {
    cat("Big Integer ('bigz') ")
    kind <- if(isM <- !is.null(nr <- attr(x, "nrow")))
      sprintf("%d x %d matrix", nr, n/nr)
    else if(n > 1) sprintf("object of length %d", n) else ""
    cat(kind,":\n", sep="")
    print(as.character(x), quote = quote, ...)
  }
  else
    cat("bigz(0)\n")
  invisible(x)
}

as.bigz <- function(a, mod = NA)
{
  if(is.null(mod)) mod <- NA
  if(inherits(a, "bigq"))
    as.bigz.bigq(a, mod)
  else
    .Call(biginteger_as, a, mod)
}

as.character.bigz <- function(x, b = 10L, ...)
{
    .Call(biginteger_as_character, x, b)
}

as.double.bigz  <- function(x,...) .Call(biginteger_as_numeric, x)
as.integer.bigz <- function(x,...) .Call(biginteger_as_integer, x)

length.bigz <- function(x) .Call(biginteger_length, x)

"length<-.bigz"<- function(x, value) .Call(biginteger_setlength, x, value)

modulus      <- function(a) UseMethod("modulus")
modulus.bigz <- function(a) attr(a, "mod")

`modulus<-`      <- function(a, value) UseMethod("modulus<-")
`modulus<-.bigz` <- function(a, value) as.bigz(a, value)


## inv <- function(a,...) UseMethod("inv")

## pow <- function(a,...) UseMethod("pow")

powm <- function(x,y, n) .Call(biginteger_powm, x,y,n)

"<.bigz"  <- function(e1, e2) .Call(biginteger_lt, e1, e2)
">.bigz"  <- function(e1, e2) .Call(biginteger_gt, e1, e2)
"<=.bigz" <- function(e1, e2) .Call(biginteger_lte, e1, e2)
">=.bigz" <- function(e1, e2) .Call(biginteger_gte, e1, e2)
"==.bigz" <- function(e1, e2) .Call(biginteger_eq, e1, e2)
"!=.bigz" <- function(e1, e2) .Call(biginteger_neq, e1, e2)

is.na.bigz <- function(x) .Call(biginteger_is_na, x)
is.finite.bigz <- function(x) !is.na.bigz(x) # otherwise all are finite
is.infinite.bigz <- function(x) rep.int(FALSE, length(x))

frexpZ <- function(x) .Call(bigI_frexp, x)

##' @title log2(Inverse of frexpZ(a))
##' @param L list(d = ., exp = .)
##' @return numeric vector
##' @author Martin Maechler
lg2.invFrexp <- function(L) {
    stopifnot(is.list(L), is.numeric(d <- L$d), is.numeric(ex <- L$exp),
	      (n <- length(d)) == length(ex))
    ex + log2(d)
}

###------------------------- 'Math' S3 group ------------------------------

## Most 'Math' group would be hard to implement --- [TODO via Rmpfr -- or stop("...via Rmpfr")?
## Fall-back: *not* implemented  {or use as.double() ??}
Math.bigz <- function(x, ...) { .NotYetImplemented() }
            ## • ‘abs’, ‘sign’, ‘sqrt’,
            ##   ‘floor’, ‘ceiling’, ‘trunc’,
            ##   ‘round’, ‘signif’
            ## • ‘exp’, ‘log’, ‘expm1’, ‘log1p’,
            ##   ‘cos’, ‘sin’, ‘tan’,
            ##   ‘acos’, ‘asin’, ‘atan’
            ##   ‘cosh’, ‘sinh’, ‘tanh’,
            ##   ‘acosh’, ‘asinh’, ‘atanh’
            ## • ‘lgamma’, ‘gamma’, ‘digamma’, ‘trigamma’
            ## • ‘cumsum’, ‘cumprod’, ‘cummax’, ‘cummin’

abs.bigz <- function(x) .Call(biginteger_abs,x)
sign.bigz <- function(x) .Call(biginteger_sgn,x)

floor.bigz <- ceiling.bigz <- function(x) x
trunc.bigz <- function(x, ...) x
round.bigz <- function(x, digits=0) {
    if(digits == 0) x
    else stop("digits != 0  is not yet implemented")
}

gamma.bigz <- function(x) factorialZ(x-1)

cumsum.bigz <- function(x) .Call(biginteger_cumsum, x)
## TODO: add cummax(), cummin(), cumprod()


log2.bigz <- function(x) .Call(biginteger_log2, x)
## not exported:
ln.bigz	 <- function(x) .Call(biginteger_log, x)
log.bigz <- function(x, base=exp(1))
{
    if(missing(base))
	ln.bigz(x)
    else
	ln.bigz(x)/log(base)
}

log10.bigz <- function(x) ln.bigz(x) / log(10)

##------------end{'Math'} group -------------------------------------


###------------------------- 'Summary' S3 group ------------------------------
##---- "max"   "min"   "range" "prod"  "sum"   "any"   "all"  -----

max.bigz <- function(...,na.rm=FALSE)
{
 .Call(biginteger_max, c.bigz(...), na.rm)
}

min.bigz <- function(...,na.rm=FALSE)
{
 .Call(biginteger_min, c.bigz(...), na.rm)
}

## range(): works automatically via  range.default() and the above min(), max()

prod.bigz <- function(..., na.rm = FALSE)
{
    X <- c.bigz(...)
   .Call(biginteger_prod, if(na.rm) X[!is.na(X)] else X)
}

sum.bigz <- function(..., na.rm = FALSE)
{
    X <- c.bigz(...)
    .Call(biginteger_sum, if(na.rm) X[!is.na(X)] else X)
}
##------------end{Summary group}------------------------------------

c.bigz <- function(..., recursive = FALSE)
{
    .Call(biginteger_c, list(...))
}

rep.bigz <- function(x,times,...)
    .Call(biginteger_rep, x, times)

# Isprime, return:
#   0 if not prime
#   1 if probably prime
#   2 if prime
isprime <- function(n,reps=40)
  {
    .Call(biginteger_is_prime, n, as.integer(reps))
  }

nextprime <- function(n) .Call(biginteger_nextprime, n)

gcdex <- function(a, b) .Call(biginteger_gcdex, a, b)

urand.bigz <- function(nb=1, size=200, seed=0)
  {
    ok <- (seed != 0)
    .Call(biginteger_rand_u,
          as.integer(nb),
          as.integer(size),
          seed,
          as.integer(ok))
  }

sizeinbase <- function(a, b=10)
{
    if(as.integer(b) < 2)
        stop("base must be >= 2")
    .Call(biginteger_sizeinbase, a, as.integer(b))
}

factorialZ  <- function(n) .Call(bigI_factorial, n)
chooseZ  <- function(n, k) .Call(bigI_choose,  n, k)

fibnum  <- function(n) .Call(bigI_fibnum,  n)
fibnum2 <- function(n) .Call(bigI_fibnum2, n)

lucnum  <- function(n) .Call(bigI_lucnum,  n)
lucnum2 <- function(n) .Call(bigI_lucnum2, n)

factorize <- function(n) .Call(factorR, as.bigz(n))

## overload as.vector
as.vector.bigz <- function(x,mode="any")
  return(x)

solve.bigz <- function(a, b,...)
  {
    if(missing(b))
      .Call(inverse_z, a)
    else
      .Call(solve_z, a, b)
  }


"[[.bigz" <- function(x, i=NA)
{
    .Call(biginteger_get_at, x, i)
}

"[[<-.bigz"<- function(dst, idx=NA, value)
{
    .Call(biginteger_set_at, dst, idx, value)
}


"[.bigz"<- function(x, i=NULL, j=NULL)
{
  .Call(matrix_get_at_z, x, i,j)
}

"[<-.bigz"<- function(dst,idx=NULL,jdx=NULL,value)
{
  .Call(matrix_set_at_z, dst, value,idx,jdx )
}


