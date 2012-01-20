library(gmp)

## want to test all these
(ops <- sapply(getGroupMembers("Ops"), getGroupMembers))

xx <- c(NaN, NA, -Inf, -123:-121, -1:2, 7:8, Inf)
(xxI <- as.bigz(xx))# Inf's and NaN's do not exist
(x <- c(NA, xx[is.finite(xx)]))
xI <- as.bigz(x)

stopifnot(identical(is.na(x), is.na(xI)),
          identical(is.finite(x), is.finite(xI)),
          identical(is.infinite(x), is.infinite(xI)))


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

opEQ <- function(OP, x, xI=as.bigz(x)) {
    if(identical(OP, `/`) || identical(OP, `%/%`) || identical(OP, `%%`)) {
        x  <- x [xn0 <- x != 0]
        xI <- xI[xn0]
    }
    EQ(mOuter(x ,, OP) -> R,
       mOuter(xI,, OP)) &
    EQ(mOuter(x,xI, OP) -> S,
       mOuter(xI,x, OP)) &
    EQ(R, S)
}

## "Compare" - works "out of the box
eqC <- lapply(sapply(ops$Compare, get),
              function(op) opEQ(op, x, xI))
stopifnot(do.call(all, eqC))

opsA <- ops$Arith

try(
    eqA <- lapply(sapply(opsA, get), function(op) opEQ(op, x, xI))
)## ^.bigz for negative input

## leave away "^" for now:
eqA <- lapply(sapply(opsA[opsA != "^"], get),
              function(op) opEQ(op, x, xI))
stopifnot(Reduce(`&`, eqA[1:3]),##  +  -  *
          eqA[["/"]])
## These are not quite ok:
eqA[["%/%"]]
eqA[["%%"]]

##-- "^" ------------
c(NA^0, 1^NA)# <- in R, these both give 1, but
N. <- as.bigz(NA)
c(N.^0, 1^N.)# both "NA" ...

## need non-negative values:
x.po0 <- x >= 0
M.pow <- opEQ(`^`, x[x.po0], xI[x.po0])# almost
M.pow[1,3] <- M.pow[2,1] <- TRUE ## <-- fix the two cases ('N.') above
stopifnot(M.pow)# all the other results match

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

xQ <- as.bigq(x)

eqC <- lapply(sapply(ops$Compare, get), function(op) opEQ(op, x, xQ))
stopifnot(Reduce(`&`, eqC))

(opsA4 <- opsA[opsA != "^" & !grepl("^%", opsA)])
eqA <- lapply(sapply(opsA4, get), function(op) opEQ(op, x, xQ))
sapply(eqA, table)
## .TRUE -.TRUE *.TRUE /.TRUE
##   100    100    100     81
##                        ^^^^ not quite
