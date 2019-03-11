/*! \file bigmod.h
 *  \brief Description of class bigmod
 *
 *  \date Created: 22/05/06
 *  \date Last modified: Time-stamp: <2019-03-10 10:30:30 (antoine)>
 *
 *  \author Immanuel Scholz
 *
 *  \note Licence: GPL
 */


#ifndef BIGMOD_HEADER_
#define BIGMOD_HEADER_ 1

#include "biginteger.h"

typedef void (*gmp_binary)(mpz_t, const mpz_t, const mpz_t);

extern "C" {
  /**
   * division
   * result = a / b
   */
  void integer_div(mpz_t result,const mpz_t a, const mpz_t b);
}



/**
 * \brief class for bigmod values. Represent any integer in Z/nZ
 *
 * Represents two biginteger: a value and a modulus. These both are used
 * to operate arithmetic functions on it. If the modulus is NA, no modulus
 * to the operation result is applied. If the value is NA, the result is always NA.
 */
class bigmod {
 private:
  /** optional source */
  biginteger * value_ptr;
  biginteger * modulus_ptr;
 
 protected:
  /** \brief  Value of our bigmod -- only references*/
  biginteger & value;
  /** \brief  modulus of our bigmod representation: value %% modulus */
  biginteger & modulus;
 
 public:

  /** keep both references value / modulus
   */
  bigmod(biginteger& value_,
	 biginteger& modulus_)  :
    value_ptr(NULL),
    modulus_ptr(NULL),
    value(value_),modulus(modulus_) {};

 /** keep references value / modulus is new object.
   */
 bigmod(biginteger& value_)  :
   value_ptr(NULL),
   modulus_ptr(new biginteger()),
   value(value_),modulus(*modulus_ptr) {};


  /**
   * create 2 new objects valus / modulus.
   */
 bigmod(const biginteger& value_,
	 const biginteger& modulus_)  :
   value_ptr(new biginteger(value_)),
   modulus_ptr(new biginteger(modulus_)),
   value(*value_ptr),modulus(*modulus_ptr) {};

 bigmod(const biginteger& value_)  :
   value_ptr(new biginteger(value_)),
   modulus_ptr(new biginteger()),
   value(*value_ptr),modulus(*modulus_ptr) {};


 bigmod()  :
   value_ptr(new biginteger()),
   modulus_ptr(new biginteger()),
   value(*value_ptr),modulus(*modulus_ptr) {};

 
  /** \brief copy operator  */
  bigmod(const bigmod & rhs) : 
    value_ptr(new biginteger(rhs.getValue())),
    modulus_ptr(new biginteger(rhs.getModulus())),
   value(*value_ptr),modulus(*modulus_ptr) {
  };


  virtual ~bigmod(){
    if(value_ptr != NULL) delete value_ptr;
    if(modulus_ptr != NULL) delete modulus_ptr;
  };

  /**
   * \brief  Return as a human readible string
   */
  std::string str(int b) const;

  /** \brief assignement operator */
  bigmod & operator= (const bigmod& rhs);

  /** \brief return sign (-1 if negative, 0 if 0; +1 if positive)
   */
  inline int sgn() const
    {
      return(mpz_sgn(getValue().getValueTemp()));
    }

  bigmod  inv () const;

 
  biginteger & getValue() {
    return value;
  }

  biginteger & getModulus() {
    return modulus;
  }

  const biginteger & getValue() const{
    return value;
  }

  const biginteger & getModulus() const {
    return modulus;
  }

};


class DefaultBigMod : public bigmod {
 private:
  /** \brief  Value of our bigmod */
  biginteger valueLocal;
  /** \brief  modulus of our bigmod representation: value %% modulus */
  biginteger modulusLocal;
 
 public:
/** \brief creator
   */
  DefaultBigMod(const biginteger& value_ = biginteger(),
	 const biginteger& modulus_ = biginteger()) :
  bigmod(valueLocal,modulusLocal),
    valueLocal(value_),modulusLocal(modulus_) {
    value = valueLocal;
    modulus = modulusLocal;
}

  /** \brief copy operator  */
 DefaultBigMod(const bigmod & rhs) :
    bigmod(valueLocal,modulusLocal),
     valueLocal(rhs.getValue()),modulusLocal(rhs.getModulus()) {
    value = valueLocal;
    modulus = modulusLocal;
}

 /** \brief copy operator  */
 DefaultBigMod(const DefaultBigMod & rhs) :
    bigmod(valueLocal,modulusLocal),
     valueLocal(rhs.getValue()),modulusLocal(rhs.getModulus()) {
    value = valueLocal;
    modulus = modulusLocal;
}
  ~DefaultBigMod(){};
 
  

};


/**
 * a bigmod that has only integer.
 */
class BigModInt : public bigmod {
 private:
   /** \brief  modulus of our bigmod representation */
  biginteger modulusLocal;
 
 public:
/** \brief creator
   */
  BigModInt(biginteger& value_) :
  bigmod(value_,modulusLocal),
    modulusLocal() {
    modulus = modulusLocal;
}

  ~BigModInt(){};
 
 
};



/** \brief comparison operator
 */
bool operator!= (const bigmod& rhs, const bigmod& lhs);

/** \brief comparison operator
 */
bool operator== (const bigmod& rhs, const bigmod& lhs);



/**
 * \brief Add two bigmods together.
 *
 * If only one has a modulus set, the result will have this
 * modulus. If both bigmods disagree with the modulus, the result will not have
 * a modulus set. If none modulus for either bigmod is set, the result will not
 * have a modulus as well.
 */
DefaultBigMod operator+(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Subtract two bigmods.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
DefaultBigMod operator-(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Multiply two bigmods.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
DefaultBigMod operator*(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Divide two bigmods   a / b  :=  a * b^(-1)
 */
DefaultBigMod div_via_inv(const bigmod& a, const bigmod& b);

/**
 * \brief Divide two bigmods.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
DefaultBigMod operator/(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Calculate the modulus (remainder) of two bigmods.
 *
 * The resulting bigmod will have set the intern modulus to
 * the value of lhs, no matter what rhs.modulus or lhs.modulus
 * was before, except if rhs and lhs has both no modulus set,
 * in which case the resulting modulus will be unset too.
 */
DefaultBigMod operator%(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief Return the power of "exp" to the base of "base" (return = base^exp).
 *
 * If both moduli are unset or unequal, this may EAT your memory alive,
 * since then the infinite "pow" is used instead of the modulus "powm".
 * You  may not try to pow a value this way with an exponent that does
 * not fit into a long value.
 *
 * For other modulus description, see operator+(bigmod, bigmod)
 */
DefaultBigMod pow(const bigmod& base, const bigmod& exp);

/**
 * \brief Return the modulo inverse to x mod m. (return = x^-1 % m)
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
DefaultBigMod inv(const bigmod& x, const bigmod& m);

/**
 * \brief Return a bigmod with value (x % m) and the intern modulus set to m.
 * Intern modulus settings of x and m are ignored.
 *
 * Do not confuse this with operator%(bigmod, bigmod).
 */
DefaultBigMod set_modulus(const bigmod& x, const bigmod& m);


biginteger get_modulus(const bigmod& b1, const bigmod& b2);
/**
 * \brief Return the greatest common divisor of both parameters
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
DefaultBigMod gcd(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief  Return the least common multiply of both parameter.
 *
 * For modulus description, see operator+(bigmod, bigmod)
 */
DefaultBigMod lcm(const bigmod& rhs, const bigmod& lhs);

/**
 * \brief function used to make any binary operation between
 * two bigmod that return a bigmod (addition substraction... )
 */
DefaultBigMod create_bigmod(const bigmod& lhs, const bigmod& rhs, gmp_binary f,
		     bool zeroRhsAllowed = true) ;

#endif
