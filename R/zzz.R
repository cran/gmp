## we "need" S4 methods for dispatch on both (x,y)  .noGenerics <- TRUE
.conflicts.OK <- TRUE

.gmpVersion <- function() .Call(R_gmp_get_version)
gmpVersion <- function()
    numeric_version(sub("^([0-9]+\\.[0-9]+\\.[0-9]+).*","\\1", .gmpVersion()))

