/*! \file bigrational.cc
 *  \brief C function for class bigrational
 *
 *  \version 1
 *
 *  \date Created: 12/12/04   
 *  \date Last modified: Time-stamp: <2005-01-10 11:21:54 lucas>
 *
 *  \author Antoine Lucas (adapted from biginteger class made by
 *                         Immanuel Scholz)
 *
 *  \note Licence: GPL
 */

#define USE_RINTERNALS
#define R_NO_REMAP   // avoid collisions with stl definitions

#include "bigrational.h"
#include <Rinternals.h>

#include <stdio.h>

using std::string;

/**
 * \brief Print value
 */
string bigrational::str() const
{
    if (numerator.isNA())
	return "NA";

    string s; // sstream seems to collide with libgmp :-(
    if (!denominator.isNA())
	s = "(";
    s += numerator.str();
    if (!denominator.isNA()) {
	s += " / ";
	s += denominator.str();
	s += ")";
    }
    return s;
}



/*
Well... I would like to have something like that...

const mpq_t& bigrational::getQValueTemp() const
{
  mpq_t val;
  mpq_init(val);
  //  mpq_t_sentry val_s(val);
  mpq_set_num(val,numerator.getValueTemp());
  if(!denominator.isNA())
    mpq_set_den(val,denominator.getValueTemp());
  return val;

}
*/

/**
 * \brief Import from double
 */
void  bigrational::setQValue(double x)
{
  mpq_t val;
  mpz_t tmpZ;

  mpq_init(val);
  mpz_init(tmpZ);

  mpq_t_sentry val_s(val);
  mpz_t_sentry tmp_s(tmpZ);

  mpq_set_d (val,x);
  mpq_get_den(tmpZ,val);
  denominator.setValue(tmpZ);

  mpq_get_num(tmpZ,val);

  numerator.setValue(tmpZ);
  /*   gmp_printf("gmp %Zd\n",den);  */

}

/*
 * \brief convert bigrational to double
 */
double  bigrational::as_double() const
{
  
  mpq_t val;
  double d;
  mpq_init( val ); 
  mpq_t_sentry val_t( val ); 
  mpq_set_num(val,numerator.getValueTemp());  
  if(!denominator.isNA())  
    mpq_set_den(val,denominator.getValueTemp());

  /*  printf("%f\n",mpq_get_d(val));*/
  d  = mpq_get_d(val);
  return(d);
}

/* \brief  simplify n/d (use of mpq_canonical)
 *
 */
void  bigrational::simplify ()
{
  
  mpq_t val;
  mpz_t den;
  mpz_t nomin;

  mpq_init( val ); 
  mpz_init(den);
  mpz_init(nomin);
      
  mpq_t_sentry val_t(val);
  mpz_t_sentry tmp_d(den);
  mpz_t_sentry tmp_n(nomin);

  mpq_set_num(val,numerator.getValueTemp());  
  if(!denominator.isNA())  
    mpq_set_den(val,denominator.getValueTemp());

  mpq_canonicalize(val);
      
  mpq_get_num(nomin ,val);
  mpq_get_den(den ,val);


}
