library(gmp)

Z1 <- as.bigz(1) ; Z1[FALSE]
Q1 <- as.bigq(1) ; Q1[FALSE]
as.bigz(0[FALSE])# failed earlier
as.bigq(0[FALSE])# failed earlier

stopifnot(identical(1L, as.integer(Z1)),
          identical(1L, as.integer(Q1)))## failed earlier

stopifnot(identical(as.bigz(1[FALSE]), Z1[FALSE]),
          identical(as.bigz(1[-1]),    Z1[-1]),
          identical(Z1[-1], rep(Z1, 0))
          , ##----------- bigq -------------
          identical(as.bigq(1[FALSE]), Q1[-1]),
          identical(Q1[FALSE], Q1[-1]),
          identical(Q1[-1], rep(Q1, 0))
          )

## FIXME:
1[0] ; Z1[0] ## the first gives empty, the second "NA"

i <- 1:9
(x <- as.bigz(i, mod = 3))
stopifnot(5*x == (5*i) %% 3)
## remove modulus "the new way" (NULL did fail):
modulus(x) <- NULL
stopifnot(identical(x, as.bigz(i %% 3)))
