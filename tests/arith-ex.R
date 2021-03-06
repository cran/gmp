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
(xxI <- as.bigz(xx))# Inf's and NaN's do not exist ==> very large integers for  +/- Inf
(x <- c(NA, xx[is.finite(xx)]))
xI <- as.bigz(x)
xQ <- as.bigq(xI)
stopifnot(identical(xI, as.bigz(xQ)),
	  identical(numerator(xQ), xI)) # numerator( <NA> )

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

## Finally (2020-06-06): mixed arithmetic works :
stopifnot(exprs = {
    isEQ(xI - xQ, c(NA, rep(0, 9)))
    isEQ(xI + xQ, 2*xI)
    isEQ(xI * xQ, x^2)
    all.equal(xQ^xI, x^x)
    ## as do mixed comparisons
    (xI == xQ)[-1]
    !(xI < xQ)[-1]
    !(xI > xQ)[-1]
    (xI >= xQ)[-1]
})

## double precision factorial() is exact up to n=22
stopifnot(factorialZ(0:22) == factorial(0:22))

## factorialZ() etc must also work when passed a bigz instead of an integer;
## till Jan.2014, they silently produced nonsense.
N <- as.bigz(n <- 3:8)
stopifnot(identical(factorialZ(N),  factorialZ(n)), factorialZ (n) == factorial(n),
          identical(chooseZ(12, N), chooseZ(12, n)), chooseZ(12,n) == choose(12,n),
          identical(fibnum (N), fibnum (n)),
          identical(fibnum2(N), fibnum2(n)),
          identical(lucnum (N), lucnum (n)),
          identical(lucnum2(N), lucnum2(n)))


## This one does *NOT* distinguish  NA and NaN -- that's wanted here
EQ1 <- function(x,y) {
    (abs(x-y) <= 1e-13*(abs(x)+abs(y)) & !(nx <- is.na(x)) & !(ny <- is.na(y))) |
        (nx & ny)
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

##' @title
##' @param OP an arithmetic OPerator such +, *,.. as R function
##' @param u  numeric vector
##' @param uI a bigz/biginteger vector,  "the same" as 'u'.
##' @return a logical  n x n matrix, say R,  R[i,j] := TRUE iff
##'   u[i] OP v[j]  are all the same when u,v vary in  {u, uI}.
##' @author Martin Maechler
opEQ <- function(OP, u, uI=as.bigz(u), eq=TRUE) {
    stopifnot(length(u) == length(uI))
    if(eq) stopifnot(isEQ(u, uI)) # should be the case when result should be all TRUE
    ##
    ## choose only some on the RHS:
    iR <-
	if(no0.R <- (identical(OP, `/`) || identical(OP, `%/%`) || identical(OP, `%%`))) {
	    ## no zero on the RHS i.e., 2nd operand
	    is.na(u) | u != 0
	} else TRUE
    ## choose only some on the LHS:
    iL <-
	if(no0.L <- (identical(OP, `^`))) {
	    ## no zero on the LHS i.e., 1st operand
	    is.na(u) | u != 0
	} else TRUE
    ##
    EQ(mOuter(u [iL],u [iR], OP) -> R,
       mOuter(uI[iL],uI[iR], OP)) &
    EQ(mOuter(u [iL],uI[iR], OP) -> S,
       mOuter(uI[iL], u[iR], OP)) &
    EQ(R, S)
}

## "Compare" - works "out of the box
eqC <- lapply(sapply(ops$Compare, get),
              function(op) opEQ(op, x, xI))
stopifnot(do.call(all, eqC))

opsA <- ops$Arith

eqA <- lapply(sapply(opsA, get), function(op) opEQ(op, x, xI))

op6 <- c("+","-", "*", "/", "%/%", "^")## << are fine - now including "^" _and_ %/% !
stopifnot(sapply(eqA, all)[op6])
## The others: now (2014-07): only %% is left: has several "wrong":
lapply(eqA[is.na(match(names(eqA), op6))], symnum)

## For example:
symnum(opEQ(`%%`, x, xI))# not all TRUE, since, e.g.,
x [3] %% x
x [3] %% xI ## (negative turned into >= 0; warning 'division by zero')

x  %% x [3]
xI %% x [3] ## (no negatives ..)


##-- "^" ------------
z1i <- 0:1
z1n <- as.double(z1i)
c(NA^0, NA^0L, z1i^NA, z1n^NA)# <- in R (<= 2011), the first and last are 1
stopifnot(isEQ(c(N.^0, N.^0L, z1i^N.), c(1,1,NA,1)),
	  isEQ(c(Nq^0, Nq^0L, z1i^Nq), c(1,1,NA,1)))

## need non-negative values:
x.po0 <- x >= 0
stopifnot(M.pow  <- opEQ(`^`, x[x.po0], xI[x.po0]))
if(FALSE)# FIXME
stopifnot(M.powQ <- opEQ(`^`, x[x.po0], xQ[x.po0]))
if(FALSE)# FIXME {z - q}
M.poIQ <- opEQ(`^`,xI[x.po0], xQ[x.po0])

## Modulo arithmetic
i <- as.bigz(-5:10, 16); i <- i[i != 0]; i
stopifnot(identical(as.integer(i), c(11:15, 1:10)))
(Ii <- 1/i )## BUG: in all versions of gmp up to 0.5-5 -- now 7 warnings pow(x, -|n|)
I2 <- i^(-1)## BUG: not considering (mod) // segmentation fault in gmp 0.5-1 {now: 7 warn..}
stopifnot(identical(Ii, I2),
	  is.na(Ii[c(2, 4, 7, 9, 11, 13, 15)]),
	  identical(Ii[c(1,3)], as.bigz(c(3,5), 16)))
(Iz <- 1/(z <- as.bigz(1:12, 13)))
stopifnot(identical(Iz, z^-1),
	  Iz == c(1, 7, 9, 10, 8, 11, 2, 5, 3, 4, 6, 12),
	  identical(modulus(Iz), as.bigz(13)))
## The first two of course give fractions:
(r1 <- as.bigz(3) / 1:12)
r2 <- as.bigz(3) / as.bigz(1:12)
stopifnot(identical(r1, r2))

## Now, the new scheme :
(iLR <- as.bigz(3, 13) / as.bigz(1:12, 13))
 ## [1] (3 %% 13)  (8 %% 13)  (1 %% 13)  (4 %% 13)  (11 %% 13) (7 %% 13)
 ## [7] (6 %% 13)  (2 %% 13)  (9 %% 13)  (12 %% 13) (5 %% 13)  (10 %% 13)
iL  <- as.bigz(3, 13) / as.bigz(1:12)
iLi <- as.bigz(3, 13) / 1:12
iR  <- as.bigz(3)     / as.bigz(1:12, 13)
iiR <-         3      / as.bigz(1:12, 13)
stopifnot(identical(iL, iLi)
          , identical(iR, iiR)
          , identical(iR, iLR)
          , identical(iL, iR)) ## failed until recently...

## whereas these two always use  divq.bigz :
(q <- as.bigz(3, 13) %/% as.bigz(1:12))
## [1] (3 %% 13) (1 %% 13) (1 %% 13) (0 %% 13) (0 %% 13) (0 %% 13)
## [7] (0 %% 13) (0 %% 13) (0 %% 13) (0 %% 13) (0 %% 13) (0 %% 13)
stopifnot(identical(q, divq.bigz(as.bigz(3, 13), 1:12)),
	  ##           ---------
	  identical(q, 3 %/% as.bigz(1:12, 13)),
	  q == c(3, 1, 1, rep(0,9)))
s <- as.bigz(3, 13) / as.bigz(1:12, 17)
## used to give
## Big Integer ('bigz') object of length 12:
##  [1] 3 1 1 0 0 0 0 0 0 0 0 0
## but now, really just  `` drops the contradicting "mod" '' ==> uses rational:
stopifnot(identical(s, r1))

##----- Z^e (modulo m) ---------------
z12 <- as.bigz(1:12,12)
stopifnot(identical(z12^1, z12), z12^0 == 1,
	  identical(z12^2, as.bigz(rep(c(1,4,9,4,1,0), 2), 12)),
	  identical(z12^3,
		    as.bigz(c(1,8,3:5,0,7:9,4,11,0), 12)),
	  identical(z12^4, z12^2),
	  identical(z12^5, z12^3),
	  identical(z12^6, z12^2),
	  identical(z12^6, (1:12) ^ as.bigz(6, 12))
	  )

for(E in 6:20) {
    ir <- as.integer(r <- z12 ^ E)
    stopifnot(identical(modulus(r), as.bigz(12)),
	      0 <= ir, ir <= 11)
}

z17 <- as.bigz(1:16, 17)
stopifnot(z17^0 == 1, identical(z17^1, z17), identical(z17^-1, iz <- 1/z17),
	  identical(z17^-2, iz^2), (iz^2) * (sq <- z17^2) == 1,
	  modulus(sq) == 17, unique(sq) == (1:8)^2 %% 17)



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
##
xQ <- as.bigq(x, 17) # == x/17 .. *are* not equal, i.e., not expecting all TRUE:
eqQ <- lapply(sapply(ops$Compare, get),
              function(op) opEQ(op, x, xQ, eq=FALSE))
lapply(eqQ, symnum)## <- symnum, for nice output

Fn <- gmp:::pow.bigq; q <- 2.3
stopifnot(inherits(e1 <- tryCatch(Fn(q,q),          error=identity), "error"),
	  inherits(e2 <- tryCatch(q ^ as.bigq(1,3), error=identity), "error"),
	  grepl("Rmpfr", e1$message),
	  identical(e1$message, e2$message))


## FIXME(2):  %% and %/%  do not work at all for  bigq
(opsA4 <- opsA[opsA != "^" & !grepl("^%", opsA)])
eqA1 <- lapply(sapply(opsA4, get), function(op) opEQ(op, x, xQ1))
sapply(eqA1, table)
## .TRUE -.TRUE *.TRUE /.TRUE
##   100    100    100     90
##                        ^^^^ (90: was 81) [not dividing by 0]

## xQ *is* different from x (apart from x[6] (and, NA x[1]))
eqA <- lapply(sapply(opsA4, get), function(op) opEQ(op, x, xQ, eq=FALSE))
lapply(eqA, symnum)


## round(x, digits) -- should work *and* be vectorized in both  (x, digits)
x1 <- as.bigq((-19:19), 10)
stopifnot(round(x1, 1) == x1)

half <- as.bigq(1, 2)
i1 <- (-19:29)
x <- half + i1
cbind(x, round(x))
rx1 <- round(x/10, 1)
stopifnot(exprs = {
    as.bigz(round(x)) %% 2 == 0
    identical(round(x) > x, i1 %% 2 == 1)
    (rx1 - x/10) * 20 == c(1,-1) # {recycling up/down}: perfect rounding to even
    (round(x/100, 2) - x/100) * 200 == c(1,-1) #  (ditto)
})
(drx1 <- asNumeric(rx1))# shows perfect round to *even*
## but double precision rounding cannot be perfect (as numbers are not exact!):
dx   <- asNumeric(x/10)
dx1  <- round(dx, 1)
dmat <- cbind(x=dx, r.x = dx1, rQx = drx1)
## shows "the picture" a bit {see Martin's vignette in CRAN package 'round'}:
noquote(cbind(apply(dmat, 2, formatC),
              ER = ifelse(abs(dx1 - drx1) > 1e-10, "*", "")))

## standard R:
rd <- round(pi*10^(-2:5), digits=7:0)
formatC(rd, digits=12, width=1)
## bigq -- show we vectorize in both x, digits
(rQ <- round(as.bigq(pi*10^(-2:5)), digits=7:0))
stopifnot(exprs = {
    as.integer(numerator  (rQ)) == 314159L
    as.integer(denominator(rQ)) == 10^(7:0)
    all.equal(asNumeric(rQ), rd, tol = 1e-15)
})


