/*! \file factor.cc
 *  \brief C function used for factorization
 *
 *  \version 1
 *
 *  \date Created: 04/12/04
 *  \date Last modified: Time-stamp: <2013-06-08 22:47:18 antoine>
 *
 *  \author Antoine Lucas (help from Immanuel Scholz) (R adaptation)
 *          Original C code from libgmp.
 *
 *  \note Licence: GPL (>= 2)
 */

#include "Rgmp.h"
#include "factorize.h"
using namespace std;


#include "factor.h"

//
// \brief function that gets values from R and send to functions
// factor
//
SEXP factorR (SEXP n)
{
  bigvec v = bigintegerR::create_bignum(n), result;
  if(v.size() > 0) {
    mpz_t val;
    mpz_init(val);
    mpz_t_sentry val_s(val);
    mpz_set(val,v[0].getValue().getValueTemp());

    int sgn = mpz_sgn(val);
    if(sgn == 0)
      error(_("Cannot factorize 0"));
    if(sgn<0)
      {
	mpz_abs(val,val);
	result.value.push_back(biginteger(-1));
      }
    //
    // function from gmplib, in demo/factorize.c
    //
    factor(val,result);
  }
  return bigintegerR::create_SEXP(result);
}
