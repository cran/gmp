#-------------------------------------------------------
#
#  Created       : 26/09/04
#  Last Modified : Time-stamp: <2004-09-27 10:58:54 lucas>
#
#  Description   : Gmp
#                  
#  Author        : Antoine Lucas
#                  antoinelucas@libertysurf.fr
#
#  Licence       : GPL 
#
#-------------------------------------------------------


#isprime <- function(n)
#  {
#    n <- as.character(n)
#    reps <- 10
#    res <- 1

#    while (res == 1)
#      {
#        res <- isprimeprob(n,reps)
        # res = 1 i.e. n probably prime
#        reps <- reps +10
#      }
#    res
#  }

# Isprime, return:
#   0 if not prime
#   1 if probably prime
#   2 if prime
isprimeprob <- function(n,reps=10)
  {
    n <- as.character(n)
    val <- 0
    .C("isprime",n,as.integer(val),as.integer(reps),
       PACKAGE= "gmp"
       )[2]
    
  }


nextprime<- function(n)
  {
    n <- as.character(n)
    val <- "0"
    .Call("nextprime",n,val,
       PACKAGE= "gmp"
       )
    
  }
gmp.add<- function(a,b)
  {
    a <- as.character(a)
    b <- as.character(b)

    .Call("gmp_add",a,b,
       PACKAGE= "gmp"
       )
    
  }
gmp.sub<- function(a,b)
  {
    a <- as.character(a)
    b <- as.character(b)

    .Call("gmp_sub",a,b,
       PACKAGE= "gmp"
       )
    
  }
gmp.mul<- function(a,b)
  {
    a <- as.character(a)
    b <- as.character(b)

    .Call("gmp_mul",a,b,
       PACKAGE= "gmp"
       )
    
  }
gmp.divq<- function(a,b)
  {
    a <- as.character(a)
    b <- as.character(b)

    .Call("gmp_fdivq",a,b,
       PACKAGE= "gmp"
       )
    
  }
gmp.divr<- function(a,b)
  {
    a <- as.character(a)
    b <- as.character(b)

    .Call("gmp_fdivr",a,b,
       PACKAGE= "gmp"
       )
    
  }
gcd<- function(a,b)
  {
    a <- as.character(a)
    b <- as.character(b)

    .Call("gmp_gcd",a,b,
       PACKAGE= "gmp"
       )
    
  }
lcm<- function(a,b)
  {
    a <- as.character(a)
    b <- as.character(b)

    .Call("gmp_lcm",a,b,
       PACKAGE= "gmp"
       )
    
  }

gmp.div <- function(a,b)
  {
    a <- as.character(a)
    b <- as.character(b)
    .C("gmp_div",a,b,as.double(0),
    PACKAGE= "gmp"
       )[3]
  }

gmp.inv <- function(a,b)
  {
    a <- as.character(a)
    b <- as.character(b)
    .Call("gmp_invert",a,b,
    PACKAGE= "gmp"
       )
  }
 
powm <- function(a,exp,modulo)
  {
    a <- as.character(a)
    exp <- as.character(exp)
    modulo <- as.character(modulo)
    .Call("gmp_powm",a,exp,modulo,
       PACKAGE= "gmp"
       )
  }

pow <- function(a,exp)
  {
    a <- as.character(a)
    .Call("gmp_pow",a,as.integer(exp),
       PACKAGE= "gmp"
       )
  }

gcdex <- function(a,b)
  {
    a <- as.character(a)
    b <- as.character(b)
    .Call("gmp_gcdex",a,b,
       PACKAGE= "gmp"
       )
  }

gmp.rand<- function(n,seed=0)
  {
    if(seed==0)
      ok <- 0
    else
      ok <- 1

    .Call("gmp_rand_u",
          as.integer(n),
          as.character(seed),
          as.integer(ok),
       PACKAGE= "gmp"
       )
  }
