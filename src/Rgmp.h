#ifndef R_gmp_HEADER
#define R_gmp_HEADER 1

// gmp.h calls cstddef with __need_size_t defined
#include <cstddef>
// and avoid the inclusion of stdlib.h
#include <cstdlib>
#include <cmath>
#include <gmp.h>



#include <R.h>
#include <Rinternals.h>
// the only thing we use from <Rdefines.h> :
#define AS_INTEGER(x) coerceVector(x,INTSXP)

#define class_P(_x_) CHAR(Rf_asChar(Rf_getAttrib(_x_, R_ClassSymbol)))

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("main", String)
#else
#define _(String) (String)
#endif
#endif
