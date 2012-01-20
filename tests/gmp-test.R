library(gmp)
## for now -- later we want *.Rout.save file --> not these details
#sessionInfo()
#packageDescription("gmp")

##
##' @title Test a *binary* function
##' @param FUN a function, such as add.bigq() ...
##' @param out string determining output class; if "str", use characters, otherwise double
##' @param x a list of "numbers"
##' @return
##' @author Antoine Lucas (& Martin Maechler)
##' @examples test(as.bigq, 0)
test <- function(FUN, out = "str", unary = FALSE,
                 x = list(23,"25", 2.3, -4, as.integer(4), 0, as.bigz(34),
                 as.bigq(32,7), as.bigz(31,45),NULL,NA,as.integer(-3)))
{
  ##    x<- list(-5,0, as.bigq(1,4),0.5,5,"10",as.integer(100), as.bigz(1000), as.bigz(3,10),NULL,NA)
  xlabs <- sapply(x, format)
  stopifnot((n <- length(x)) >= 1)
  if(out == "str") {
      sortie <- as.character
      res <- ""
      error <- "error"
  } else {
      sortie <- as.double
      res <- 0
      error <- NA
  }
  nr <- if(unary) 1 else n
  res <- matrix(res, nr, n)

  for(i in 1:nr)
    for(j in 1:n) {
      e <- if(unary) try(FUN(x[[j]]),        silent = TRUE) else
                     try(FUN(x[[i]],x[[j]]), silent = TRUE)
      if(inherits(e, "try-error"))
        e <- error
      else if(length(e) == 0)
        e <- numeric()

      res[i,j] <- sortie(e)[1]
    }
  res <- as.data.frame(res)
  ##    dimnames(res)[[1]]<- labels
  names(res) <- xlabs
  res
}## end{test}


allfunctionid <- c("as.bigz","add.bigz","sub.bigz","mul.bigz",
                   "divq.bigz","div.bigz","mod.bigz","pow.bigz",
                   "inv.bigz","gcd.bigz","lcm.bigz",
                   "add.bigq","sub.bigq","div.bigq","as.bigq","gcdex",
		   "max.bigq","max.bigz","min.bigq","min.bigz")
unaryfunctionid <- c("log.bigz","log2.bigz","log10.bigz","c.bigz",
                     "isprime","nextprime",
                   "sizeinbase","fibnum","fibnum2","lucnum","lucnum2",
                     "factorize","abs")

options(width = 125)

gmp.NS <- asNamespace("gmp")# also get namespace *hidden* functions, i.e. methods:
for(fid in allfunctionid)
  {
    cat ("------------------------------------------\n", fid, "\n\n", sep="")
    FUN <- get(fid, envir=gmp.NS, mode="function")
    print(res <- test(FUN))
  }

##==============================================================================

for(fid in unaryfunctionid)
  {
    cat ("------------------------------------------\n", fid, "\n\n", sep="")
    FUN <- get(fid, envir = gmp.NS, mode="function")
    print(res <- test(FUN, unary=TRUE))
  }

##==============================================================================

###----------- matrix -----------------------------
x  <- matrix(1:6,3)
stopifnot(identical(as.bigz(x), matrix(as.bigz(as.vector(x)), 3)),
          dim(x) == 3:2,
          dim(x) == dim(ym   <- as.bigz(x, 6:1)),
          dim(x) == dim(ymr  <- as.bigz(x, 4:6)),
          dim(x) == dim(ymc  <- as.bigz(x, 4)),
          dim(x) == dim(ymq  <- as.bigq(x)),
          dim(x) == dim(y    <- as.bigq(x, 6:1))
          ,
          apply(ym,1,max) == 1:3,
          apply(ym,2,min) == c(1,1,0))

x %*% t(x)

ym %*% t(ym)
ym %*% t(ymr)
ymc %*% t(ymc)
ymq %*% t(ymq)
y %*% t(y)

dd <- dim(D  <- diag(1:4))
stopifnot(dd == dim(Dmq  <- as.bigq(D)),
          dd == dim(Dm   <- as.bigz(D,6:1)),
          dd == dim(Dmr  <- as.bigz(D,7)),
          dd == dim(Dmc  <- as.bigz(D,4)),
          TRUE)
solve(D)
solve(Dmq)
solve(Dmr)
try(solve(Dmc))# Error: argument has no inverse
try(solve(Dm)) # Error: System is singular

(D.D <- try(D %*% t(Dm)))# now [>= Jan.2012] works too

