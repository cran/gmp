library(gmp)

## for reference  (==> *not* using *.Rout.save here!)
sessionInfo()
packageDescription("gmp")

##' an (x == y) which gives TRUE also when both are NA:
isEQ <- function(x,y)  (x == y) | (is.na(x) & is.na(y))

## want to test all these
(ops <- sapply(getGroupMembers("Ops"), getGroupMembers))

N. <- as.bigz(NA)
Nq <- as.bigq(NA)
stopifnot(identical(Nq, as.bigq(N.)),
	  identical(N., as.bigz(Nq)))# used to fail

xx <- c(NaN, NA, -Inf, -123:-121, -1:2, 7:8, Inf)
(xxI <- as.bigz(xx))# Inf's and NaN's do not exist
(x <- c(NA, xx[is.finite(xx)]))
xI <- as.bigz(x)
xQ <- as.bigq(xI)
stopifnot(identical(xI, as.bigz(xQ)))

stopifnot(isEQ(x, as.integer(x)), isEQ(x, xI), isEQ(x, xQ),
	  identical(xQ, as.bigq(x)),
	  identical(is.na(x), is.na(xI)), identical(is.na(x), is.na(xQ)),
	  identical(is.finite(x), is.finite(xI)),
	  identical(is.finite(x), is.finite(xQ)),
	  identical(is.infinite(x), is.infinite(xI)),
	  identical(is.infinite(x), is.infinite(xQ)),
## The next 4 all failed till 2012-05-05:
	  isEQ(x, as.integer(xI)),
	  isEQ(x, as.integer(xQ)),
	  isEQ(x, as.numeric(xI)),
	  isEQ(x, as.numeric(xQ)),
          TRUE)

## double precision factorial() is exact up to n=22
stopifnot(factorialZ(0:22) == factorial(0:22))

## This one does *NOT* distinguish  NA and NaN -- that's wanted here
EQ1 <- function(x,y) {
    (abs(x-y) <= 1e-13*(abs(x)+abs(y)) & !(nx <- is.na(x)) & !(ny <- is.na(y))) | (nx & ny)
}
stopifnot(EQ1(x, xI))
EQ <- function(x,y) mapply(EQ1, x, y, USE.NAMES=FALSE)

## a version of outer() that should work also with these objects
mOuter <- function(X, Y=X, FUN, ...) {
    lapply(seq_along(X), function(i) FUN(X[i], Y, ...))
}

matOuter <- function(X, Y=X, FUN, ...) {
    t(array(unlist(mOuter(X, Y, FUN, ...)),
	    dim = c(length(Y), length(X))))
}


opEQ <- function(OP, x, xI=as.bigz(x)) {
    iR <-
	if(no0.R <- (identical(OP, `/`) || identical(OP, `%/%`) || identical(OP, `%%`))) {
	    ## no zero on the RHS i.e., 2nd operand
	    is.na(x) | x != 0
	} else TRUE
    iL <-
	if(no0.L <- (identical(OP, `^`))) {
	    ## no zero on the LHS i.e., 1st operand
	    is.na(x) | x != 0
	} else TRUE

    EQ(mOuter(x [iL],x [iR], OP) -> R,
       mOuter(xI[iL],xI[iR], OP)) &
    EQ(mOuter(x [iL],xI[iR], OP) -> S,
       mOuter(xI[iL], x[iR], OP)) &
    EQ(R, S)
}

## "Compare" - works "out of the box
eqC <- lapply(sapply(ops$Compare, get),
              function(op) opEQ(op, x, xI))
stopifnot(do.call(all, eqC))

opsA <- ops$Arith

eqA <- lapply(sapply(opsA, get), function(op) opEQ(op, x, xI))


op5 <- c("+","-", "*","/", "^") ## These are fine - now including "^" !
stopifnot(sapply(eqA, all)[op5])
## The others:  %/% and %% have several "wrong"
lapply(eqA[is.na(match(names(eqA), op5))], symnum)

##-- "^" ------------
z1i <- 0:1
z1n <- as.double(z1i)
c(NA^0, NA^0L, z1i^NA, z1n^NA)# <- in R (<= 2011), the first and last are 1
stopifnot(isEQ(c(N.^0, N.^0L, z1i^N.), c(1,1,NA,1)),
	  isEQ(c(Nq^0, Nq^0L, z1i^Nq), c(1,1,NA,1)))

## need non-negative values:
x.po0 <- x >= 0
stopifnot(M.pow  <- opEQ(`^`, x[x.po0], xI[x.po0]))
stopifnot(M.powQ <- opEQ(`^`, x[x.po0], xQ[x.po0]))
if(FALSE)# FIXME {z - q}
M.poIQ <- opEQ(`^`,xI[x.po0], xQ[x.po0])


##--- Log()s -------------------------
(ex <- c(outer(c(2,5,10), 10^(1:3))))# 20 .. 10'000
stopifnot(dim(L <- outer(as.bigz(2:4), ex, `^`)) == c(3, length(ex)))
l2 <- array(log2(L), dim = dim(L))
lnL <- log(L)
a.EQ <- function(x,y, tol=1e-15, ...) all.equal(x,y, tol=tol, ...)
stopifnot(a.EQ(l2[1,], ex),
          a.EQ(l2[3,], 2*ex),
          a.EQ(log(L, 8), lnL/log(8)),
          a.EQ(c(l2), lnL/log(2)))


###------------------ bigq --------------------------------


xQ1 <- as.bigq(x, 1)
eqC <- lapply(sapply(ops$Compare, get), function(op) opEQ(op, x, xQ1))
stopifnot(Reduce(`&`, eqC))##
## but this fails:
xQ <- as.bigq(x, 17)
eqC <- lapply(sapply(ops$Compare, get), function(op) opEQ(op, x, xQ))
lapply(eqC, symnum)
sapply(eqC, all) # not one ...  {FIXME}?

(opsA4 <- opsA[opsA != "^" & !grepl("^%", opsA)])

eqA1 <- lapply(sapply(opsA4, get), function(op) opEQ(op, x, xQ1))
sapply(eqA1, table)
## .TRUE -.TRUE *.TRUE /.TRUE
##   100    100    100     90
##                        ^^^^ (90: was 81) .. not quite
eqA <- lapply(sapply(opsA4, get), function(op) opEQ(op, x, xQ))
sapply(eqA, table)
##        +  -  *  /
## FALSE 80 80 64 64
## TRUE  20 20 36 26 ('26' was '17') <<--- oy oy ... way to go!

## and  %%  and %/%  do not work at all for  bigq  {but do for numeric!} --- FIXME
