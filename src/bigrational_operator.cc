/*! \file bigrational_operator.cc
 *  \brief C operator for bigrational
 *
 *  \version 1
 *
 *  \date Created: 12/12/04   
 *  \date Last modified: Time-stamp: <2005-02-06 19:22:13 antoine>
 *
 *  \author Antoine Lucas (adapted from biginteger class made by
 *                         Immanuel Scholz)
 *
 *  \note Licence: GPL
 */


#define USE_RINTERNALS
#define R_NO_REMAP // avoid collisions with stl definitions


#include "bigrational.h"

#include <Rinternals.h>

#include "stdio.h"

namespace
{
  /*
    // Use this to clear mpz_t structs at end-of-function automatically
    struct mpz_t_sentry {
	mpz_t& value;
	mpz_t_sentry(mpz_t& v): value(v) {}
	~mpz_t_sentry() {mpz_clear(value);}
    };
    // Use this to clear mpq_t structs at end-of-function automatically
    struct mpq_t_sentry {
	mpq_t& value;
	mpq_t_sentry(mpq_t& v): value(v) {}
	~mpq_t_sentry() {mpq_clear(value);}
    };
  */
    typedef void (*gmpq_binary)(mpq_t, const mpq_t, const mpq_t);
    // Create a bigrational from a binary combination of two other bigrationals
    bigrational create_bigrational(const bigrational& lhs, const bigrational& rhs, gmpq_binary f,  bool zeroRhsAllowed = true) {
	if (lhs.numerator.isNA() || rhs.numerator.isNA())
	    return bigrational();

	if (!zeroRhsAllowed && (rhs.as_double() == 0))
	  {
	    printf("rhs %f\n",rhs.as_double());
	    Rf_error("division by zero");
	  }

	mpq_t val;
	mpq_init(val);
	mpq_t_sentry val_s(val);
	mpq_set_num(val,lhs.numerator.getValueTemp());
	if(!lhs.denominator.isNA())
	  mpq_set_den(val,lhs.denominator.getValueTemp());


	mpq_t val2;
	mpq_init(val2);
	mpq_t_sentry val2_s(val2);
	mpq_set_num(val2,rhs.numerator.getValueTemp());
	if(!rhs.denominator.isNA())
	  mpq_set_den(val2,rhs.denominator.getValueTemp());


	mpz_t den;
	mpz_init(den);
	mpz_t_sentry den_s(den);

	mpz_t nomin;
	mpz_init(nomin);
	mpz_t_sentry nomin_s(nomin);

        f(val,val,val2);

	/* Simplify numerator and denominator */	
	mpq_canonicalize(val);

	/*	gmp_printf("gmp %Qd\n",val);*/

	mpq_get_num(nomin ,val);
	mpq_get_den(den ,val);
	
	return bigrational(nomin, den);
    }
}


/**
 * \brief Return  a + b
 */
bigrational operator+(const bigrational& lhs, const bigrational& rhs)
{

    return create_bigrational(lhs, rhs, mpq_add);
}

/**
 * \brief Return  a - b
 */
bigrational operator-(const bigrational& lhs, const bigrational& rhs)
{
    return create_bigrational(lhs, rhs, mpq_sub);
}

/**
 * \brief Return  a * b
 */
bigrational operator*(const bigrational& lhs, const bigrational& rhs)
{
    return create_bigrational(lhs, rhs, mpq_mul);
}

/**
 * \brief Return  a / b
 */
bigrational operator/(const bigrational& lhs, const bigrational& rhs)
{
    return create_bigrational(lhs, rhs, mpq_div, false);
}


/**
 * \brief Well... an heritage from biginteger class, this should be
 * integrated earlier... put denominator & simplify if there is not.
 */
bigrational set_denominator(const bigrational& n, const bigrational& d)
{

  if(d.numerator.as_double() == 0)
    {
      printf("Division by 0 \n");
      return bigrational();
    }
  
  /* All this to use mpq_canonicalize i.e: simplify n/d */ 
  mpq_t valnomin;
  mpq_t valden;
  mpz_t den;
  mpz_t nomin;
  
  mpq_init(valnomin);
  mpq_init(valden);
  mpz_init(den);
  mpz_init(nomin);
      
  mpq_t_sentry val1_s(valnomin);
  mpq_t_sentry val2_s(valden);
  mpz_t_sentry tmp_d(den);
  mpz_t_sentry tmp_n(nomin);
  if(d.numerator.isNA())
    printf("fuck:!!\n");
  mpq_set_num(valnomin,n.numerator.getValueTemp());
  mpq_set_num(valden,d.numerator.getValueTemp());
  if(!n.denominator.isNA())
    mpq_set_den(valnomin,n.denominator.getValueTemp());
  if(!d.denominator.isNA())
    mpq_set_den(valden,d.denominator.getValueTemp());
  /* simplify n/d */
  
  mpq_div(valnomin ,valnomin,valden);
  //      mpq_canonicalize(valnomin);
  
  mpq_get_num(nomin ,valnomin);
  mpq_get_den(den ,valnomin);

  return bigrational(nomin,den);
}





