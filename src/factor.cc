/*! \file factor.cc
 *  \brief C function used for factorization
 *
 *  \version 1
 *
 *  \date Created: 04/12/04   
 *  \date Last modified: Time-stamp: <2004-12-04 10:13:41 antoine>
 *
 *  \author Antoine Lucas (help from Immanuel Scholz) (R adaptation)
 *          Original C code from libgmp.
 *
 *  \note Licence: GPL
 */


#define USE_RINTERNALS

#define R_NO_REMAP    	// avoid collisions with stl definitions

#include <math.h>
#include <gmp.h>

#include <R.h>
#include <Rdefines.h>

#undef PROTECT
#undef UNPROTECT
#define PROTECT(x) Rf_protect(x)	// but use some handy defs anyways
#define UNPROTECT(x) Rf_unprotect(x)
#undef coerceVector
#define coerceVector             Rf_coerceVector

#include "biginteger.h"

#include <stdio.h>

#include <vector>
#include <algorithm>
using namespace std;
#include "bigintegerR.h"


extern "C"
{

  SEXP factorR (SEXP n);

}


int flag_verbose = 0;

static unsigned add[] = {4, 2, 4, 2, 4, 6, 2, 6};

/** brief Use this to clear mpz_t structs at end-of-function automatically
 */
struct mpz_t_sentry {
  mpz_t& value;
  mpz_t_sentry(mpz_t& v): value(v) {}
  ~mpz_t_sentry() {mpz_clear(value);}
};


void factor_using_division (mpz_t t, unsigned int limit,  vector<bigmod> * result)
{
  mpz_t q, r;
  unsigned long int f;
  int ai;
  unsigned *addv = add;
  unsigned int failures;

  if (flag_verbose)
    {
      printf ("[trial division (%u)] ", limit);
      fflush (stdout);
    }

  mpz_init (q);
  mpz_init (r);

  f = mpz_scan1 (t, 0);
  mpz_div_2exp (t, t, f);
  while (f)
    {

      (*result).push_back(bigmod());
      (*result)[(*result).size()-1].value.setValue(2);
      //      printf ("2 ");
      fflush (stdout);
      --f;
    }

  for (;;)
    {
      mpz_tdiv_qr_ui (q, r, t, 3);
      if (mpz_cmp_ui (r, 0) != 0)
        break;
      mpz_set (t, q);
      //printf ("3 ");
      fflush (stdout);
      (*result).push_back(bigmod());
      (*result)[(*result).size()-1].value.setValue(3);
     }

  for (;;)
    {
      mpz_tdiv_qr_ui (q, r, t, 5);
      if (mpz_cmp_ui (r, 0) != 0)
        break;
      mpz_set (t, q);
      //printf ("5 ");
      fflush (stdout);
      (*result).push_back(bigmod());
      (*result)[(*result).size()-1].value.setValue(5);
    }

  failures = 0;
  f = 7;
  ai = 0;
  while (mpz_cmp_ui (t, 1) != 0)
    {
      mpz_tdiv_qr_ui (q, r, t, f);
      if (mpz_cmp_ui (r, 0) != 0)
        {
          f += addv[ai];
          if (mpz_cmp_ui (q, f) < 0)
            break;
          ai = (ai + 1) & 7;
          failures++;
          if (failures > limit)
            break;
        }
      else
        {
          mpz_swap (t, q);
          //printf ("%lu ", f);
          fflush (stdout);
          failures = 0;
          (*result).push_back(bigmod());
          (*result)[(*result).size()-1].value.setValue(f);
        }
    }

  mpz_clear (q);
  mpz_clear (r);
}

void
factor_using_division_2kp (mpz_t t, unsigned int limit, unsigned long p,  vector<bigmod> * result)
{
  mpz_t r;
  mpz_t f;
  unsigned int k;

  mpz_init (r);
  mpz_init_set_ui (f, 2 * p);
  mpz_add_ui (f, f, 1);
  for (k = 1; k < limit; k++)
    {
      mpz_tdiv_r (r, t, f);
      while (mpz_cmp_ui (r, 0) == 0)
        {
          mpz_tdiv_q (t, t, f);
          mpz_tdiv_r (r, t, f);
          (*result).push_back(bigmod());
          (*result)[(*result).size()-1].value.setValue(f);
          //mpz_out_str (stdout, 10, f);
          fflush (stdout);
          fputc (' ', stdout);
        }
      mpz_add_ui (f, f, 2 * p);
    }

  mpz_clear (f);
  mpz_clear (r);
}

void
factor_using_pollard_rho (mpz_t n, int a_int, unsigned long p, vector<bigmod> * result)
{
  mpz_t x, x1, y, P;
  mpz_t a;
  mpz_t g;
  mpz_t t1, t2;
  int k, l, c, i;

  if (flag_verbose)
    {
      printf ("[pollard-rho (%d)] ", a_int);
      fflush (stdout);
    }

  mpz_init (g);
  mpz_init (t1);
  mpz_init (t2);

  mpz_init_set_si (a, a_int);
  mpz_init_set_si (y, 2);
  mpz_init_set_si (x, 2);
  mpz_init_set_si (x1, 2);
  k = 1;
  l = 1;
  mpz_init_set_ui (P, 1);
  c = 0;

  while (mpz_cmp_ui (n, 1) != 0)
    {
S2:
      if (p != 0)
        {
          mpz_powm_ui (x, x, p, n); mpz_add (x, x, a);
        }
      else
        {
          mpz_mul (x, x, x); mpz_add (x, x, a); mpz_mod (x, x, n);
        }
      mpz_sub (t1, x1, x); mpz_mul (t2, P, t1); mpz_mod (P, t2, n);
      c++;
      if (c == 20)
        {
          c = 0;
          mpz_gcd (g, P, n);
          if (mpz_cmp_ui (g, 1) != 0)
            goto S4;
          mpz_set (y, x);
        }
S3:
      k--;
      if (k > 0)
        goto S2;

      mpz_gcd (g, P, n);
      if (mpz_cmp_ui (g, 1) != 0)
        goto S4;

      mpz_set (x1, x);
      k = l;
      l = 2 * l;
      for (i = 0; i < k; i++)
        {
          if (p != 0)
            {
              mpz_powm_ui (x, x, p, n); mpz_add (x, x, a);
            }
          else
            {
              mpz_mul (x, x, x); mpz_add (x, x, a); mpz_mod (x, x, n);
            }
        }
      mpz_set (y, x);
      c = 0;
      goto S2;
S4:
      do
        {
          if (p != 0)
            {
              mpz_powm_ui (y, y, p, n); mpz_add (y, y, a); 
            }
          else
            {
              mpz_mul (y, y, y); mpz_add (y, y, a); mpz_mod (y, y, n);
            }
          mpz_sub (t1, x1, y); mpz_gcd (g, t1, n);
        }
      while (mpz_cmp_ui (g, 1) == 0);

      if (!mpz_probab_prime_p (g, 3))
        {
          do
            {
              mp_limb_t a_limb;
              mpn_random (&a_limb, (mp_size_t) 1);
              a_int = (int) a_limb;
            }
          while (a_int == -2 || a_int == 0);

          if (flag_verbose)
            {
              printf ("[composite factor--restarting pollard-rho] ");
              fflush (stdout);
            }
          factor_using_pollard_rho (g, a_int, p,result);
          break;
        }
      else
        {
          //mpz_out_str (stdout, 10, g);
          (*result).push_back(biginteger());
          (*result)[(*result).size()-1].value.setValue(g);

          fflush (stdout);
          fputc (' ', stdout);
        }
      mpz_div (n, n, g);
      mpz_mod (x, x, n);
      mpz_mod (x1, x1, n);
      mpz_mod (y, y, n);
      if (mpz_probab_prime_p (n, 3))
        {
          (*result).push_back(biginteger());
          (*result)[(*result).size()-1].value.setValue(n);
          //mpz_out_str (stdout, 10, n);
          fflush (stdout);
          fputc (' ', stdout);
          break;
        }
    }

  mpz_clear (g);
  mpz_clear (P);
  mpz_clear (t2);
  mpz_clear (t1);
  mpz_clear (a);
  mpz_clear (x1);
  mpz_clear (x);
  mpz_clear (y);
}

void
factor (mpz_t t, unsigned long p,  vector<bigmod> * result)
{
  unsigned int division_limit;

  if (mpz_sgn (t) == 0)
    return;

  /* Set the trial division limit according the size of t.  */
  division_limit = mpz_sizeinbase (t, 2);
  if (division_limit > 1000)
    division_limit = 1000 * 1000;
  else
    division_limit = division_limit * division_limit;

  if (p != 0)
    factor_using_division_2kp (t, division_limit / 10, p,result);
  else
    factor_using_division (t, division_limit,result);

  if (mpz_cmp_ui (t, 1) != 0)
    {
      if (flag_verbose)
        {
          printf ("[is number prime?] ");
          fflush (stdout);
        }
      if (mpz_probab_prime_p (t, 3))
        {
          (*result).push_back(biginteger());
          (*result)[(*result).size()-1].value.setValue(t);
          //mpz_out_str (stdout, 10, t);

        }
      else
        factor_using_pollard_rho (t, 1, p,result);
    }
}

/*
 * \brief function that gets values from R and send to functions
 * factor
 */
SEXP factorR (SEXP n)
{
  vector<bigmod> v = bigintegerR::create_bignum(n);
  vector<bigmod> result;
    
  result.reserve(0);
  
  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
  mpz_set(val,v[0].value.getValueTemp());
  factor(val,0,&result);


  return bigintegerR::create_SEXP(result);

}
