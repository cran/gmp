
#include <gmp.h>
#include <R.h>
#include <Rdefines.h>

#include <string>

/**
 * A big integer. Actually a wrapper for mpz_t to work with plus 
 * some special stuff.
 *
 * The biginteger special state "NA" means, no value is assigned. 
 * This does not mean, the internal state is not constructed, but 
 * the value explicit is "not available".
 */


#include "biginteger.h"

class bigrational
{
    
public:
  /**
   * The rational = n/d .
   */

  biginteger numerator;
  biginteger denominator;

  /**
   * Construct a bigrational from 2 bigintegers 
   */
  bigrational(const biginteger& numerator_ = biginteger(), 
	 const biginteger& denominator_ = biginteger()) : numerator(numerator_),denominator(denominator_) {}



  //  const mpq_t& getQValueTemp() const ;
  /*
   * \brief Set a value from double
   */
  void setQValue(double x);

  /* \brief  simplify n/d (use of mpq_canonical)
   *
   */
  void  simplify ();

  /*
   * \brief Get a double value
   */
  double as_double() const;

  /* Macro to get a mpq (mpr) from a biginteger (bigr) 
   * Example:
   *   mpq_t val;
   *   GET_MPQ ( br, val,val_tmp)
   */
  #define GET_MPQ( bigr,mpr,tmp) mpq_init( mpr ); mpq_t_sentry tmp( mpr ); mpq_set_num(mpr,bigr.numerator.getValueTemp());  if(!bigr.denominator.isNA())  mpq_set_den(mpr,bigr.denominator.getValueTemp());
  
    
  /**
   * Return as a human readible string
   */
  std::string str(int b) const;
};

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



bigrational set_denominator(const bigrational& x, const bigrational& m);


/**
 * Add two bigrationals together.
 *
 * If only one has a modulus set, the result will have this
 * modulus. If both bigrationals disagree with the modulus, the result will not have
 * a modulus set. If none modulus for either bigrational is set, the result will not
 * have a modulus as well.
 */
bigrational operator+(const bigrational& rhs, const bigrational& lhs);

/**
 * Subtract two bigrationals.
 *
 * For modulus description, see operator+(bigrational, bigrational)
 */
bigrational operator-(const bigrational& rhs, const bigrational& lhs);

/**
 * Multiply two bigrationals.
 *
 * For modulus description, see operator+(bigrational, bigrational)
 */
bigrational operator*(const bigrational& rhs, const bigrational& lhs);

/**
 * Divide two bigrationals.
 *
 * For modulus description, see operator+(bigrational, bigrational)
 */
bigrational operator/(const bigrational& rhs, const bigrational& lhs);

/**
 * Calculate the modulus (remainder) of two bigrationals.
 *
 * The resulting bigrational will have set the intern modulus to 
 * the value of lhs, no matter what rhs.modulus or lhs.modulus
 * was before, except if rhs and lhs has both no modulus set,
 * in which case the resulting modulus will be unset too.
 */
//bigrational operator%(const bigrational& rhs, const bigrational& lhs);

/**
 * Return the power of "exp" to the base of "base" (return = base^exp).
 *
 * If both moduli are unset or unequal, this may EAT your memory alive, 
 * since then the infinite "pow" is used instead of the modulus "powm". 
 * You  may not try to pow a value this way with an exponent that does
 * not fit into a long value.
 *
 * For other modulus description, see operator+(bigrational, bigrational)
 */
//bigrational pow(const bigrational& base, const bigrational& exp);

/**
 * Return the modulo inverse to x mod m. (return = x^-1 % m)
 *
 * For modulus description, see operator+(bigrational, bigrational)
 */
//bigrational inv(const bigrational& x, const bigrational& m);


