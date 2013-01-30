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


add.bigq <- `+.bigq` <- function(e1, e2) .Call(bigrational_add, e1, e2)

sub.bigq<- `-.bigq` <- function(e1, e2=NULL) {
  if(is.null(e2))
    .Call(bigrational_sub,0,e1)
  else
    .Call(bigrational_sub,e1,e2)
}

mul.bigq <- `*.bigq` <- function(e1, e2) .Call(bigrational_mul, e1, e2)

"/.bigq" <- div.bigq <- function(e1, e2) .Call(bigrational_div, e1, e2)

## note: Here, the 2nd argument must be big integer, not rational:
"^.bigq" <- pow.bigq <- function(e1, e2) .Call(bigrational_pow, e1, as.bigz(e2))

print.bigq <- function(x, quote = FALSE, initLine = TRUE, ...)
{
  if((n <- length(x)) > 0) {
    if(initLine) {
      cat("Big Rational ('bigq') ")
      kind <- if(isM <- !is.null(nr <- attr(x, "nrow")))
        sprintf("%d x %d matrix", nr, n/nr)
      else if(n > 1) sprintf("object of length %d", n) else ""
      cat(kind,":\n", sep="")
    }
    print(as.character(x), quote = quote, ...)
  }
  else
    cat("bigq(0)\n")
  invisible(x)
}

as.bigq <- function(n, d=1)
{
    .Call(bigrational_as, n, d)
}

as.character.bigq <- function(x, b = 10L, ...)
{
    .Call(bigrational_as_character, x, b)
}

formatN.bigq	<- function(x, ...) {
    r <- as.character(x, ...)
    if(any(iI <- is.whole.bigq(x)))
	r[iI] <- paste0(r[iI],"/1")
    r
}

as.double.bigq<- function(x,...) .Call(bigrational_as_numeric, x)
## maybe sub-optimal, but at least "R-consistent" in warnings/errors...:
as.integer.bigq<- function(x,...) as.integer(.Call(bigrational_as_numeric, x))

.bigq2num <- function(x) {
    r <- .Call(bigrational_as_numeric, x)
    if(!is.null(d <- dim(x))) dim(r) <- d
    r
}
setMethod("asNumeric", "bigq", .bigq2num)


denominator <- function(x) {
  r <- .Call(bigrational_den,x)
  if(!is.null(d <- dim(x))) dim(r) <- d
  r
}

"denominator<-" <- function(x,value)
  as.bigq(numerator(x),value)

numerator <- function(x) {
    r <- .Call(bigrational_num,x)
    if(!is.null(d <- dim(x))) dim(r) <- d
    r
}

"numerator<-" <- function(x,value)
  as.bigq(value,denominator(x))


as.bigz.bigq <- function(a, mod = NA)
{
  ## "FIXME":  considerably faster in C / C++
  if(any(ina <- is.na.bigq(a))) {
    r <- as.bigz(rep.int(NA, length(a)))
    if(any(ii <- !ina)) {
	a <- a[ii]
	r[ii] <- as.bigz(numerator(a) %/% denominator(a), mod[ii])
    }
    attr(r,"nrow") <- attr(a, "nrow")
    r
  }
  else # no NA's
    as.bigz(numerator(a) %/% denominator(a), mod)
}

length.bigq<- function(x) .Call(bigrational_length, x)
`length<-.bigq` <- function(x, value) .Call(bigrational_setlength, x, value)


"<.bigq"  <- function(e1,e2) .Call(bigrational_lt,  e1, e2)
">.bigq"  <- function(e1,e2) .Call(bigrational_gt,  e1, e2)
"<=.bigq" <- function(e1,e2) .Call(bigrational_lte, e1, e2)
">=.bigq" <- function(e1,e2) .Call(bigrational_gte, e1, e2)
"==.bigq" <- function(e1,e2) .Call(bigrational_eq,  e1, e2)
"!=.bigq" <- function(e1,e2) .Call(bigrational_neq, e1, e2)

is.na.bigq <- function(x) .Call(bigrational_is_na, x)
is.whole.bigq <- function(x) .Call(bigrational_is_int, x)
is.finite.bigq <- function(x) !is.na.bigq(x) # otherwise all are finite
is.infinite.bigq <- function(x) rep.int(FALSE, length(x))

###  <bigz> o <bigq>  --- really dispatch on two arguments --> use S4
if(FALSE) { ## not working really --- see also ./matrix-prods.R
##                cannot use as.bigz() yet which is only defined in ./bigz.R
setOldClass("bigz")#, prototype=as.bigz(integer()))
setOldClass("bigq")#, prototype=as.bigq(integer()))

setMethod("Ops", signature(e1 = "bigq", e2 = "bigz"),
	  function(e1, e2) callGeneric(e1, as.bigq(e2)))
setMethod("Ops", signature(e1 = "bigz", e2 = "bigq"),
	  function(e1, e2) callGeneric(as.bigq(e1), e2))
}

###------------------------- 'Math' S3 group ------------------------------

## Most 'Math' group would be hard to implement --- [TODO via Rmpfr -- or stop("...via Rmpfr")?
## Fall-back: *not* implemented  {or use as.double() ??}
Math.bigq <- function(x, ...) { .NotYetImplemented() }


abs.bigq <- function(x) {
    numerator(x) <- abs(numerator(x))
    x
}

sign.bigq <- function(x) sign(numerator(x))

## TODO
trunc.bigq <- function(x, ...) .NotYetImplemented()
floor.bigq   <- function(x) .NotYetImplemented()
ceiling.bigq <- function(x) .NotYetImplemented()
round.bigq <- function(x, digits = 0) {
    .NotYetImplemented()
}

cumsum.bigq <- function(x) .Call(bigrational_cumsum, x)
## TODO: add cummax(), cummin(), cumprod()

## FIXME: implement  log() etc --- see ./biginteger.R

##------------end{'Math'} group -------------------------------------




c.bigq <- function(..., recursive = FALSE) {
    .Call(bigrational_c, list(...))
}

rep.bigq <- function(x,times,...) .Call(bigrational_rep,x,times)

solve.bigq <- function(a,b,...)
  {
    if(missing(b))
      .Call(inverse_q,a)
    else
      .Call(solve_q,a,b)
  }


`[[.bigq`<- function(x, i=NA)
{
    .Call(bigrational_get_at, x, i)
}

`[[<-.bigq` <- function(x, i=NA, value)
{
    .Call(bigrational_set_at, x, i, value)
}


`[.bigq` <- function(x,i=NULL,j=NULL, drop=TRUE)
{
  .Call(matrix_get_at_q, x, i,j)
}


`[<-.bigq` <- function(x,i=NULL,j=NULL,value)
{
  .Call(matrix_set_at_q, x, value,i,j )
}


max.bigq <- function(...,na.rm=FALSE)
{
 .Call(bigrational_max, c.bigq(...), na.rm)
}

min.bigq <- function(...,na.rm=FALSE)
{
 .Call(bigrational_min, c.bigq(...), na.rm)
}

sum.bigq <- function(..., na.rm = FALSE)
{
    X <- c.bigq(...)
   .Call(bigrational_sum, if(na.rm) X[!is.na(X)] else X)
}

prod.bigq <- function(..., na.rm = FALSE)
{
    X <- c.bigq(...)
   .Call(bigrational_prod, if(na.rm) X[!is.na(X)] else X)
}



