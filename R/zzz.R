.noGenerics <- TRUE
.conflicts.OK <- TRUE

.onLoad <- .First.lib <- function(lib, pkg)
{
    library.dynam("gmp", pkg, lib)

}

