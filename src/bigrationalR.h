/*! \file bigrationalR.h
 *  \brief header for C++ functions that deals with both bigrational and R
 *
 *  \version 1
 *
 *  \date Created: 2005   
 *  \date Last modified: Time-stamp: <2006-05-26 14:42:00 antoine>
 *
 *
 *  \note Licence: GPL
 */
#ifndef BIG_RATIONALRRRR_
#define BIG_RATIONALRRRR_ 1

#include "bigvec_q.h"

typedef bigrational (*bigrational_binary_fn)(const bigrational&, const bigrational&);
typedef bool (*bigrational_logical_binary_fn)(const bigrational&, const bigrational&);

extern "C"
{
  /**
   * \brief Addition of a and b
   */
  SEXP bigrational_add(SEXP a, SEXP b);

  /**
   * \brief Subtraction of a and b
   */
  SEXP bigrational_sub(SEXP a, SEXP b);

  /**
   * \brief Multiplication of a and b
   */
  SEXP bigrational_mul(SEXP a, SEXP b);

  /**
   * \brief Quotient of a / b
   */
  SEXP bigrational_div(SEXP a, SEXP b);

  /**
   * \brief Return Numerator of a 
   */
  SEXP bigrational_num(SEXP a);

  /**
   * \brief Return Denominator of a 
   */
  SEXP bigrational_den(SEXP a);
   
  /**
   * \brief Return from vector a all elements specified in vector b
   */
  SEXP bigrational_get_at(SEXP a, SEXP b);
    
  /**
   * \brief Return a vector with the values from src specified by 
   * idx to sequentiell values from "value".
   */
  SEXP bigrational_set_at(SEXP src, SEXP idx, SEXP value);
    
  /**
   * \brief Convert from one or 2 long value or a string or
   * bigz into bigrational.
   *
   */
  SEXP bigrational_as(SEXP n, SEXP d);

  /**
   * \brief Convert from a bigrational vector to a character string vector.
   */
  SEXP bigrational_as_character(SEXP a, SEXP b);

  /**
   * \brief Convert from a bigrational vector to a real vector.
   */
  SEXP bigrational_as_numeric(SEXP a);

  /**
   * \brief Return the length of the vector
   */
  SEXP bigrational_length(SEXP a);

  /**
   * \brief Returns a resized vector cut at end or filled with NA.
   */
  SEXP bigrational_setlength(SEXP vec, SEXP value);

  /**
   * \brief Return whether the parameter is NA
   */
  SEXP bigrational_is_na(SEXP a);

  /**
   * \brief Return whether a < b
   */
  SEXP bigrational_lt(SEXP a, SEXP b);
    
  /**
   * \brief Return whether a > b
   */
  SEXP bigrational_gt(SEXP a, SEXP b);
    
  /**
   * \brief Return whether a <= b
   */
  SEXP bigrational_lte(SEXP a, SEXP b);
    
  /**
   * \brief Return whether a >= b
   */
  SEXP bigrational_gte(SEXP a, SEXP b);
    
  /**
   * \brief Return whether a == b
   */
  SEXP bigrational_eq(SEXP a, SEXP b);
    
  /**
   * \brief Return whether a != b
   */
  SEXP bigrational_neq(SEXP a, SEXP b);

  /**
   * \brief For function c()
   */
  SEXP bigrational_c(SEXP args) ;

  /** \brief for function cbind
   */ 
  SEXP bigrational_cbind(SEXP args) ;
  /**
   * \brief Create vector as n times x
   */
  SEXP bigrational_rep(SEXP x, SEXP times) ;

  /**
   * \brief Return max
   */
  SEXP bigrational_max(SEXP a, SEXP narm);

  /**
   * \brief Return min
   */
  SEXP bigrational_min(SEXP a, SEXP narm);

  /**
   * \brief Return cumsum
   */
  SEXP bigrational_cumsum(SEXP a);
  /**
   * \brief Return prod
   */
  SEXP bigrational_prod(SEXP a);


}


/**
 * \brief set of function usefull for manipulation of SEXP and bigvec_q
 *
 */
namespace bigrationalR{

  bigvec_q create_vector(SEXP param);

  bigvec_q create_bignum(SEXP param);

  SEXP create_SEXP(const bigvec_q & v);
  
  SEXP bigrational_binary_operation(SEXP a, SEXP b, bigrational_binary_fn f);

  SEXP bigrational_logical_binary_operation(SEXP a, SEXP b, bigrational_logical_binary_fn f);

  typedef void (*gmpq_binary)(mpq_t, const mpq_t, const mpq_t);

  bigrational create_bigrational(const bigrational& lhs, const bigrational& rhs, gmpq_binary f,  bool zeroRhsAllowed = true);

  bool lt(const bigrational& lhs, const bigrational& rhs);

  bool gt(const bigrational& lhs, const bigrational& rhs) ;

  bool lte(const bigrational& lhs, const bigrational& rhs) ;

  bool gte(const bigrational& lhs, const bigrational& rhs) ;

  bool eq(const bigrational& lhs, const bigrational& rhs);

  bool neq(const bigrational& lhs, const bigrational& rhs);


}


#endif
