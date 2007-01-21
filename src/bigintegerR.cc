/************************************************************/
/*! \file bigintegerR.cc
 *  \brief C function to interface R and libgmp with big integer values
 *
 *  \version 1
 *
 *  \date Created: 27/10/04   
 *  \date Last modified: Time-stamp: <2006-06-14 22:59:27 antoine>
 *
 *  \author Immanuel Scholz (help from A. Lucas)
 *
 *  \note Licence: GPL
 */

#define USE_RINTERNALS

#define R_NO_REMAP 			// avoid collisions with stl definitions

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


#include "bigintegerR.h"
#include "matrix.h"

#include <stdio.h>

#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;


/* Globals variables */

static gmp_randstate_t seed_state;
static int seed_init=0;

namespace bigintegerR
{
  // \brief create a vector of bigvecs, all without a modulus.
  bigvec create_vector(SEXP param) {
    switch (TYPEOF(param)) {
    case RAWSXP:
      {	
	// deserialise the vector. first int is the size.
	bigvec v;
	const char* raw = (char*)RAW(param);
	int pos = sizeof(int); // position in raw[]. Starting after header.
	int sizevec = ((int*)raw)[0];
	//std::cout << "nb element a lire " << sizevec << std::endl;
	v.value.resize(sizevec);
	for (int i = 0; i < sizevec; ++i) {
	  v.value[i] = biginteger(&raw[pos]);
	  pos += v.value[i].raw_size(); // increment number of bytes read.
	}
	return v;
      }
    case REALSXP:
      {
	double* d = REAL(param);
	//bigvec v(d,d+LENGTH(param));
	bigvec v;	
	v.value.resize(LENGTH(param));
	for (int j = 0; j < LENGTH(param); ++j)
	  v.value[j] = d[j];
	return v;
      }
    case INTSXP:
    case LGLSXP:
      {
	int* i = INTEGER(param);
	//bigvec v(i,i+LENGTH(param));
	bigvec v;
	v.value.resize(LENGTH(param));
		
	for (int j = 0; j < LENGTH(param); ++j)
	    v.value[j] = i[j];

	return v;
      }
    case STRSXP:
      {
	bigvec v;
	v.value.resize(LENGTH(param));
	for (int i = 0; i < LENGTH(param); ++i) {
	  if (STRING_ELT(param,i) == NA_STRING)
	    v.value[i]= biginteger();
	  else
	    v.value[i]=biginteger(std::string(CHAR(STRING_ELT(param,i))));
	}
	return v;
      }
    default:
      return bigvec();
    }
  }

  bigvec create_bignum(SEXP param) {
    SEXP modName;
    PROTECT(modName = Rf_allocVector(STRSXP,1));
    SET_STRING_ELT(modName, 0, Rf_mkChar("mod"));
    SEXP modAttr = Rf_getAttrib(param, modName);
    UNPROTECT(1);

    SEXP dimAttr, dimName;
    PROTECT(dimAttr);
    PROTECT(dimName = Rf_allocVector(STRSXP,1));    
    SET_STRING_ELT(dimName, 0, Rf_mkChar("nrow"));
    dimAttr = Rf_getAttrib(param, dimName);
    UNPROTECT(2);

    // try to catch biz-nrow dimension value
    //std::cout << "import value" << std::endl;
    bigvec v = bigintegerR::create_vector(param);

    if (TYPEOF(dimAttr) == INTSXP) 
      v.nrow = INTEGER(dimAttr)[0];
    else
      {
	// catch to get std matrix dimensions value
	PROTECT(dimAttr);
	PROTECT(dimName = Rf_allocVector(STRSXP,1));    
	SET_STRING_ELT(dimName, 0, Rf_mkChar("dim"));
	dimAttr = Rf_getAttrib(param, dimName);
	UNPROTECT(2);
	if (TYPEOF(dimAttr) == INTSXP) 
	  v.nrow = INTEGER(dimAttr)[0];
	else
	  v.nrow = 0;
      }

    if (TYPEOF(modAttr) != NILSXP) 
      {
	//std::cout << "import value" << std::endl;
	v.modulus = bigintegerR::create_vector(modAttr).value;
      }
    return v;
	
  }

  std::vector<int> create_int(SEXP param) {
    switch (TYPEOF(param)) {
    case REALSXP:
      {
	double* d = REAL(param);
	// copy vector manually to avoid stupid conversation warning in STL-code :-/
	vector<int> v;
	v.reserve(LENGTH(param));
	for (int i = 0; i < LENGTH(param); ++i)
	  v.push_back(static_cast<int>(d[i]));
	return v;
	//return vector<int>(d, d+LENGTH(param));
      }
    case INTSXP:
    case LGLSXP:
      {
	int* i = INTEGER(param);

	return std::vector<int>(i, i+LENGTH(param));
      }
    default:
      return std::vector<int>();
    }
  }


  SEXP create_SEXP(const std::vector<biginteger>& v) 
  {
    SEXP ans;
    unsigned int i;
    int size = sizeof(int); // starting with vector-size-header
    for (i = 0; i < v.size(); ++i)
      size += v[i].raw_size(); // adding each bigint's needed size
    PROTECT(ans = Rf_allocVector(RAWSXP, size ));
    char* r = (char*)(RAW(ans));
    ((int*)(r))[0] = v.size(); // first int is vector-size-header
    int pos = sizeof(int); // current position in r[] (starting after vector-size-header)
    for (i = 0; i < v.size(); ++i)
      pos += v[i].as_raw(&r[pos]);
    UNPROTECT(1);
    return(ans);
  }

  SEXP create_SEXP(const bigvec& v) {

    SEXP ans;

    PROTECT(ans = create_SEXP(v.value));
	
    // set the class attribute to "bigz"
    SEXP rexpName;
    PROTECT(rexpName = Rf_allocVector(STRSXP, 1));
    SET_STRING_ELT(rexpName, 0, Rf_mkChar("bigz"));
    Rf_setAttrib(ans, R_ClassSymbol, rexpName);

    // set the dim attribute to "bigz"
    
    if(v.nrow != 0 )
      {
	//	    SEXP nrowName;
	PROTECT(rexpName = Rf_allocVector(STRSXP, 1));
	SET_STRING_ELT(rexpName, 0, Rf_mkChar("nrow"));
	SEXP nRow;
	PROTECT(nRow = Rf_allocVector(INTSXP, 1));
	INTEGER(nRow)[0] = (int) v.nrow;
	Rf_setAttrib(ans, rexpName,nRow);
	UNPROTECT(2);
      }
    
    // set the mod attribute
    if (v.modulus.size()>0)
      {		
	SEXP modAttr = create_SEXP(v.modulus); 
	PROTECT(rexpName = Rf_allocVector(STRSXP,1));
	SET_STRING_ELT(rexpName, 0, Rf_mkChar("mod"));
	Rf_setAttrib(ans, rexpName, modAttr);
	UNPROTECT(1);
      }
    UNPROTECT(2);

    return ans;
  }

  /**
   * \brief Main function of doing a binary operation on bigintegers.
   * It calls a function argument for doing the correct thing.
   * This could also be written as a class functor (template)
   * to save one function call, but then code bloat will happen.
   */
  SEXP biginteger_binary_operation(SEXP a, SEXP b, biginteger_binary_fn f)
  {
    bigvec va = bigintegerR::create_bignum(a);
    bigvec vb = bigintegerR::create_bignum(b);
    if (va.value.empty() || vb.value.empty())
      Rf_error("argument must not be an empty list.");
    bigvec result;
    int size = std::max(va.value.size(), vb.value.size());
    result.value.reserve(size);
    for (int i = 0; i < size; ++i)
	result.push_back(f(va[i%va.size()], vb[i%vb.size()]));

    result.nrow = matrixz::checkDims(va.nrow,vb.nrow) ;
    return bigintegerR::create_SEXP(result);
  }


  SEXP biginteger_logical_binary_operation(SEXP a, SEXP b, biginteger_logical_binary_fn f)
  {
    bigvec va = bigintegerR::create_bignum(a);
    bigvec vb = bigintegerR::create_bignum(b);
    if (va.value.empty() || vb.value.empty())
      Rf_error("argument must not be an empty list.");
    int size = max(va.value.size(), vb.value.size());
    //	int sizemod = max(va.modulus.size(), vb.modulus.size());
    SEXP ans;
    PROTECT(ans = Rf_allocVector(LGLSXP, size));
    /* TODO: this kind of situation 5 == (5 %% 17)*/
    for (int i = 0; i < size; ++i) {
      biginteger am = va.value[i % va.value.size()];
      biginteger bm = vb.value[i % vb.value.size()];
      if (am.isNA() || bm.isNA())
	LOGICAL(ans)[i] = NA_LOGICAL;
      else
	LOGICAL(ans)[i] = f(am, bm) ? 1 : 0;
    }


    int nrow = matrixz::checkDims(va.nrow,vb.nrow) ;

    // Add dimension parameter when available
    if(nrow>0)
      {
	SEXP dimName,dimVal;
	PROTECT(dimName = Rf_allocVector(STRSXP, 1));
	PROTECT(dimVal = Rf_allocVector(INTSXP, 2));
	SET_STRING_ELT(dimName, 0, Rf_mkChar("dim"));
	INTEGER(dimVal)[0] = (int) nrow;
	INTEGER(dimVal)[1] = (int) size / nrow;
	Rf_setAttrib(ans, dimName,dimVal);
	UNPROTECT(2);
      }

    UNPROTECT(1);
    return ans;
  }

  bool lt(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) < 0;}
  bool gt(const biginteger& lhs, const biginteger& rhs) 
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) > 0;}
  bool lte(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) <= 0;}
  bool gte(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) >= 0;}
  bool eq(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) == 0;}
  bool neq(const biginteger& lhs, const biginteger& rhs)
  {return mpz_cmp(lhs.getValueTemp(), rhs.getValueTemp()) != 0;}

}

/* End of namespace bigintegerR*/

SEXP biginteger_add (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator+);}
SEXP biginteger_sub (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator-);}
SEXP biginteger_mul (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator*);}
SEXP biginteger_div (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator/);}
SEXP biginteger_mod (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,operator%);}
SEXP biginteger_pow (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,pow);}
SEXP biginteger_inv (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,inv);}
SEXP biginteger_gcd (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,gcd);}
SEXP biginteger_lcm (SEXP a, SEXP b) {return bigintegerR::biginteger_binary_operation(a,b,lcm);}
SEXP biginteger_as (SEXP a, SEXP mod) {return bigintegerR::biginteger_binary_operation(a,mod,set_modulus);}

SEXP biginteger_lt (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::lt);}
SEXP biginteger_gt (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::gt);}
SEXP biginteger_lte (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::lte);}
SEXP biginteger_gte (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::gte);}
SEXP biginteger_eq (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::eq);}
SEXP biginteger_neq (SEXP a, SEXP b) {return bigintegerR::biginteger_logical_binary_operation(a,b,bigintegerR::neq);}



SEXP biginteger_as_character(SEXP a, SEXP b)
{
  bigvec v = bigintegerR::create_bignum(a);
  SEXP ans;
  int base;
  base = *INTEGER(b);
  if (base<2 || base>36)
    Rf_error("select a base between 2 and 36");

  PROTECT(ans = Rf_allocVector(STRSXP, v.size()));
  for (unsigned int i = 0; i < v.size(); ++i)
    SET_STRING_ELT(ans, i, Rf_mkChar(v.str(i,base).c_str()));
  // matrix part
  if(v.nrow > 0)
    {
      SEXP rexpName,nRow;
      PROTECT(rexpName = Rf_allocVector(STRSXP, 1));
      SET_STRING_ELT(rexpName, 0, Rf_mkChar("dim"));
      PROTECT(nRow = Rf_allocVector(INTSXP, 2));
      INTEGER(nRow)[0] = v.nrow;
      INTEGER(nRow)[1] = v.value.size() / v.nrow;
      Rf_setAttrib(ans, rexpName,nRow);
      UNPROTECT(2);

    }

  UNPROTECT(1);
  return ans;
}

SEXP biginteger_as_numeric(SEXP a)
{
  bigvec v = bigintegerR::create_bignum(a);
  SEXP ans;
  PROTECT(ans = Rf_allocVector(REALSXP,v.size()));
  for (unsigned int i = 0; i < v.size(); ++i)
    REAL(ans)[i] = v.value[i].as_double();
  UNPROTECT(1);
  return ans;
}

SEXP biginteger_get_at(SEXP a, SEXP b)
{
  //result = a [b]

  bigvec va = bigintegerR::create_bignum(a);
  return(bigintegerR::create_SEXP(bigintegerR::biginteger_get_at_C(va,b)));

}

bigvec bigintegerR::biginteger_get_at_C(bigvec va,SEXP b)
{
  vector<int> vb = bigintegerR::create_int(b);
  bigvec result;
  // logical: b = true/false
  if (TYPEOF(b) == LGLSXP) {
    for (unsigned int i = 0; i < va.size(); ++i)
      if (vb[i%vb.size()])
	{
	  //std::cout << "cas LOGIC "<< std::endl;
	  result.push_back(va[i]);
	}
    return result;
  } 
  else {
    vb.erase(remove(vb.begin(), vb.end(), 0), vb.end()); // remove all zeroes
    if (vb.empty())
      return bigvec();

    // case: a[-b]
    if (vb[0] < 0) {
      //std::cout << "cas ngatif" << std::cout;
      for (vector<int>::iterator it = vb.begin(); it != vb.end(); ++it)
	if (*it > 0)
	  Rf_error("only 0's may mix with negative subscripts");
	else if (-(*it)-1 >= (int)va.size())
	  Rf_error("subscript out of bounds");

      // TODO: This is optimized for large va.size and small vb.size.
      // Maybe add a condition to use a different approach for large vb's
      result.value.reserve(va.size()-vb.size());
      for (int i = 0; i < (int)va.size(); ++i)
	if (find(vb.begin(), vb.end(), -i-1) == vb.end())
	  {
	    result.push_back(va[i]);
	  }
    }
    else {
      // standard case: a[b] with b: integers
      result.value.reserve(vb.size());
      for (vector<int>::iterator it = vb.begin(); it != vb.end(); ++it) {
	if (*it <= 0)
	  Rf_error("only 0's may mix with negative subscripts");
	if (*it <= (int)va.size())
	  {
	    //std::cout << "on sort " << va.value[(*it)-1].str(10) << std::endl;
	    result.push_back(va[(*it)-1]);
	  }
	else
	  result.push_back(bigmod()); // NA for out of range's
      }
    }
  }
  return (result);
}

SEXP biginteger_set_at(SEXP src, SEXP idx, SEXP value)
{

  // return = ( src[idx] <- value )
  bigvec result = bigintegerR::create_bignum(src);
  bigvec vvalue = bigintegerR::create_bignum(value);
  vector<int> vidx = bigintegerR::create_int(idx);

  //case: logicals
  if (TYPEOF(idx) == LGLSXP) {
    int pos = 0;
    for (unsigned int i = 0; i < result.size(); ++i)
      if (vidx[i%vidx.size()])
	result.set(i, vvalue[pos++%vvalue.size()]);
    return bigintegerR::create_SEXP(result);
  }
  else {
    vidx.erase(remove(vidx.begin(), vidx.end(), 0), vidx.end()); // remove all zeroes
    if (vidx.empty())
      return bigintegerR::create_SEXP(result);
    // return = (src[-idx] <- value)
    if (vidx[0] < 0) {
      for (vector<int>::iterator it = vidx.begin(); it != vidx.end(); ++it)
	if (*it > 0)
	  Rf_error("only 0's may mix with negative subscripts");
	else if (-(*it)-1 >= (int)result.size())
	  Rf_error("subscript out of bounds");
      int pos = 0;
      for (int i = 0; i < (int)result.size(); ++i)
	if (find(vidx.begin(), vidx.end(), -i-1) == vidx.end())
	  result.set(i, vvalue[pos++%vvalue.size()]);
    }
    //standard case: return = (src[idx] <- value) with idx: positive integer
    else {
      // finding maximum to resize vector if needed
      int maximum = INT_MIN;
      for (vector<int>::iterator it = vidx.begin(); it != vidx.end(); ++it)
	maximum = max(maximum, *it);
      if (maximum > (int)result.size())
	result.resize(maximum);
      int pos = 0;
      for (vector<int>::iterator it = vidx.begin(); it != vidx.end(); ++it) {
	if (*it < 0)
	  Rf_error("only 0's may mix with negative subscripts");
	result.set((*it)-1,vvalue[pos++%vvalue.size()]);
      }
    }
  }
  return bigintegerR::create_SEXP(result);
}

SEXP biginteger_length(SEXP a)
{
  SEXP ans;
  PROTECT(ans = Rf_allocVector(INTSXP,1));
  INTEGER(ans)[0] = bigintegerR::create_bignum(a).size();
  UNPROTECT(1);
  return ans;
}

SEXP biginteger_setlength(SEXP vec, SEXP value)
{
  int len = 0;
  switch (TYPEOF(value)) {
  case INTSXP:
  case LGLSXP:
    if (LENGTH(value) != 1)
      Rf_error("invalid second argument");
    len = *INTEGER(value);
    if (len < 0)
      Rf_error("vector size cannot be negative");
    else if (len == NA_INTEGER)
      Rf_error("vector size cannot be NA");
    break;
  case REALSXP:
    if (LENGTH(value) != 1)
      Rf_error("invalid second argument");
    len = (int)*REAL(value);
    if (len < 0)
      Rf_error("vector size cannot be negative");
    else if (! (R_FINITE (len ) ))
      Rf_error("vector size cannot be NA, NaN of Inf");
    break;
  case STRSXP:
    // dunno why R spits out this strange error on "length(foo) <- -1"
    // but I always follow the holy standard ;-)
    Rf_error("negative length vectors are not allowed");
  default:
    Rf_error("invalid second argument");
  }
  bigvec v =bigintegerR::create_bignum(vec);
  v.resize(len);
  return bigintegerR::create_SEXP(v);
}

SEXP biginteger_is_na(SEXP a) 
{
  bigvec v = bigintegerR::create_bignum(a);
  SEXP ans;
  PROTECT(ans = Rf_allocVector(LGLSXP, v.size()));
  for (unsigned int i = 0; i < v.size(); ++i)
    LOGICAL(ans)[i] = v[i].value.isNA();
  UNPROTECT(1);
  return ans;
}


SEXP biginteger_sgn(SEXP a) 
{
  bigvec v = bigintegerR::create_bignum(a);
  SEXP ans;
  PROTECT(ans = Rf_allocVector(INTSXP, v.size()));
  for (unsigned int i = 0; i < v.size(); ++i)
    INTEGER(ans)[i] = mpz_sgn(v[i].value.getValueTemp());
  UNPROTECT(1);
  return ans;
}

SEXP biginteger_c(SEXP args) 
{
  //if(TYPEOF( args ) != LISTSXP)
  //  Rf_error("should be a list");

  int i=0,j=0; 
  bigvec result;
  bigvec v;

  for(i =0; i<LENGTH(args);i++)
    {
      v = bigintegerR::create_bignum(VECTOR_ELT(args,i));
      for(j=0; j< (int)v.size() ; j++)
	result.push_back(v[j]);
      v.clear();
    }

  return bigintegerR::create_SEXP(result);
}

SEXP biginteger_cbind(SEXP args) 
{
  int i=0,j=0,nrow=0; 
  bigvec result;
  bigvec v;

  result = bigintegerR::create_bignum(VECTOR_ELT(args,0));
  if(result.nrow ==0)
    result.nrow == result.size();

  for(i =1; i<LENGTH(args);i++)
    {
      v = bigintegerR::create_bignum(VECTOR_ELT(args,i));
      for(j=0; j< (int)v.size() ; j++)
	result.push_back(v[j]);
      v.clear();
    }

  return bigintegerR::create_SEXP(result);
}


SEXP biginteger_rep(SEXP x, SEXP times) 
{
  bigvec v = bigintegerR::create_bignum(x);
  bigvec result;
  int i,j,rep;
  
  PROTECT(times= AS_INTEGER(times));
  rep = INTEGER(times)[0];
  UNPROTECT(1);

  result.value.reserve(v.size()*rep);
  for(i = 0 ; i< rep ; i++)
    for(j = 0 ; j < (int)v.size() ; j++)
      result.push_back(v[j]);
  
  return bigintegerR::create_SEXP(result);
}



SEXP biginteger_is_prime(SEXP a, SEXP reps) 
{
  bigvec v = bigintegerR::create_bignum(a);
  vector<int> vb = bigintegerR::create_int(reps);
  unsigned int i;
  SEXP ans;
  PROTECT(ans = Rf_allocVector(INTSXP, v.size()));
  if(v.size() == vb.size())
    for (i = 0; i < v.size(); ++i)
      INTEGER(ans)[i] = v[i].value.isprime(vb[i]);
  else
    for (i = 0; i < v.size(); ++i)
      INTEGER(ans)[i] = v[i].value.isprime(vb[0]);
  UNPROTECT(1);
  return ans;
}
 

SEXP biginteger_nextprime(SEXP a) 
{
  bigvec v =bigintegerR::create_bignum(a);
  bigvec result;
    
  result.value.reserve(v.size());
    
  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);

  for (unsigned int i = 0; i < v.size(); ++i)
    {
      mpz_nextprime(val,v[i].value.getValueTemp());
      result.push_back(bigmod(val));
    }
    
  return bigintegerR::create_SEXP(result);
    
}

SEXP biginteger_abs(SEXP a) 
{
  bigvec v =bigintegerR::create_bignum(a);
  bigvec result;
    
  result.value.reserve(v.size());
    
  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);

  for (unsigned int i = 0; i < v.size(); ++i)
    {
      mpz_abs(val,v[i].value.getValueTemp());
      result.push_back(bigmod());
      result[i].value.setValue(val);
    }
    
  return bigintegerR::create_SEXP(result);
    
}


/** @brief Bezoult coefficients: compute g,s and t as as + bt = g 
 *  @param a BigInteger
 *  @param b BigInteger
 */
SEXP biginteger_gcdex(SEXP a, SEXP b) 
{
  bigvec va = bigintegerR::create_bignum(a);
  bigvec vb = bigintegerR::create_bignum(b);
  bigvec result;
  unsigned int i;

  if(va.size() != vb.size())
    return bigintegerR::create_SEXP(bigvec());
    

  result.value.reserve(3*va.size());

  mpz_t g;
  mpz_t s;
  mpz_t t;
  mpz_init(g);
  mpz_init(s);
  mpz_init(t);
  mpz_t_sentry val_g(g);
  mpz_t_sentry val_s(s);
  mpz_t_sentry val_t(t);

  for(i=0; i<va.size();i++)
    {
      mpz_gcdext (g,s,t,va[i].value.getValueTemp(),vb[i].value.getValueTemp());
      result.value.push_back(biginteger(g)); // Hem... not very elegant !
      result.value.push_back(biginteger(s));
      result.value.push_back(biginteger(t));
      /*      result[i*3].value.setValue(g);
      result[i*3+1].value.setValue(s);
      result[i*3+2].value.setValue(t);*/
	
    }
  return bigintegerR::create_SEXP(result);

}

/** @brief Random number generation
    \note If seed is not initialised: generation of a new seed
    @param nb  Integer: number of number to generate
    @param length Integer number will be of length 2^length
    @param newseed Integer, seed initialisation (if exists)
    @param ok Integer 1: seed generation 0 not
*/

SEXP biginteger_rand_u (SEXP nb ,SEXP length,SEXP newseed, SEXP ok)
{
  mpz_t    bz;
  int i,flag,len,size;
  bigvec result;
   
  //extern int seed_init;
  //extern gmp_randstate_t seed_state;


  /* store input data into appropriate mode */
  bigvec va = bigintegerR::create_bignum(newseed);
  PROTECT (ok = AS_INTEGER(ok));
  PROTECT (length = AS_INTEGER(length));
  PROTECT (nb = AS_INTEGER(nb));
  flag = INTEGER(ok)[0];
  len = INTEGER(length)[0];
  size = INTEGER(nb)[0];
  UNPROTECT(3);

  result.value.reserve(size);
  
  /* Random seed initialisation */

  if(seed_init==0)
    {
      gmp_randinit_default(seed_state);
      printf("Seed default initialisation\n");
    }
  if(flag == 1)
    {
      gmp_randseed(seed_state,va[0].value.getValueTemp());
      printf("Seed initialisation\n");
    }

  seed_init = 1;

  mpz_init (bz);
  mpz_t_sentry val_s(bz);
  
  for(i= 0; i<size; i++)
    {
      /*  Random number generation  */
      mpz_urandomb(bz,seed_state,len);
      result.push_back(bigmod(bz));
    }
  return bigintegerR::create_SEXP(result);
}


/** @brief biginteger_sizeinbase return 
 *  @param x BigInteger
 *  @param base BigInteger
 */
SEXP biginteger_sizeinbase(SEXP x, SEXP base) 
{
  bigvec vx = bigintegerR::create_bignum(x);

  int i,basesize;
  size_t taille;

  SEXP ans;
    
  PROTECT (base = AS_INTEGER(base));
  basesize=INTEGER(base)[0];
  UNPROTECT(1);

  PROTECT(ans = Rf_allocVector(INTSXP,vx.size()));

  for(i=0; i<(int)vx.size();i++)
    {
      taille = mpz_sizeinbase(vx[i].value.getValueTemp(),basesize);
      INTEGER(ans)[i]=taille;
      //	printf("%d \n" ,taille);
    }
  UNPROTECT(1);

  return ans;

}


/** @brief fibnum return nth Fibonacci number 
 *  @param n integer
 */
SEXP fibnum(SEXP n) 
{
  bigvec result;
  int num;

  if(LENGTH(n) > 0)
    {
      PROTECT (n = AS_INTEGER(n));
      num = INTEGER(n)[0];
      UNPROTECT(1);

      if(num<0)
	Rf_error("argument must be positive");

      mpz_t val;
      mpz_init(val);
      mpz_t_sentry val_s(val);
  
      mpz_fib_ui(val,num);
      result.push_back(bigmod(val));
      //      result[0].value.setValue(val);
    }
  else
    Rf_error("argument must not be an empty list");

  return bigintegerR::create_SEXP(result);

}

/** @brief fibnum2 return nth and n-1th Fibonacci number 
 *  @param n integer
 */
SEXP fibnum2(SEXP n) 
{
  bigvec result;
  int num;

  if(LENGTH(n) > 0)
    {
      PROTECT (n = AS_INTEGER(n));
      num = INTEGER(n)[0];
      UNPROTECT(1);

      if(num<0)
	Rf_error("argument must be positive");

      result.value.reserve(1);

      mpz_t val;
      mpz_init(val);
      mpz_t_sentry val_s(val);
      mpz_t val2;
      mpz_init(val2);
      mpz_t_sentry val_s2(val2);
      
      mpz_fib2_ui(val,val2,num);
      result.push_back(bigmod(val2));
      result.push_back(bigmod(val));

    }
  else
    Rf_error("argument must not be an empty list");

  return bigintegerR::create_SEXP(result);

}

/** @brief lucnum return nth lucas number 
 *  @param n integer
 */
SEXP lucnum(SEXP n) 
{
  bigvec result;
  int num;

  if(LENGTH(n) > 0)
    {
      PROTECT (n = AS_INTEGER(n));
      num = INTEGER(n)[0];
      UNPROTECT(1);

      if(num<0)
	Rf_error("argument must be positive");

      mpz_t val;
      mpz_init(val);
      mpz_t_sentry val_s(val);
  
      mpz_lucnum_ui(val,num);
      result.push_back(bigmod(val));
    }
  else
    Rf_error("argument must not be an empty list");

  return bigintegerR::create_SEXP(result);

}

/** @brief lucnum2 return nth and n-1th lucas number 
 *  @param n integer
 */
SEXP lucnum2(SEXP n) 
{
  bigvec result;
  int num;

  if(LENGTH(n) > 0)
    {
      PROTECT (n = AS_INTEGER(n));
      num = INTEGER(n)[0];
      UNPROTECT(1);

      if(num<0)
	Rf_error("argument must be positive");

      mpz_t val;
      mpz_init(val);
      mpz_t_sentry val_s(val);
      mpz_t val2;
      mpz_init(val2);
      mpz_t_sentry val_s2(val2);
      
      mpz_lucnum2_ui(val,val2,num);
      result.push_back(bigmod(val2));
      result.push_back(bigmod(val));
    }
  else
    Rf_error("argument must not be an empty list");

  return bigintegerR::create_SEXP(result);

}



//Return max

SEXP biginteger_max(SEXP a, SEXP narm)
{

  bigvec result;

  bigvec va = bigintegerR::create_bignum(a);

  if( ! va.size())
    return bigintegerR::create_SEXP(result);

  unsigned int maximum = 0;

  PROTECT (a = AS_INTEGER(a));
  int na_remove = INTEGER(a)[0];
  UNPROTECT(1);


  for(unsigned int i = 1 ; i < va.size(); ++i)
    {
      if(va.value[i].isNA() && 
	 ( ! na_remove) )
	return(bigintegerR::create_SEXP(result));
      else
	if(!(va.value[i] <  va.value[maximum] ))
	  maximum = i; // if va.value[maximum = 0] is NA => false for the "<" => maximum changed = good
    }

  result.push_back(va.value[maximum]);


  // now the modulus !
  if(va.modulus.size() == 1)
    result.modulus.push_back(va.modulus[0]);

  if(va.modulus.size()>1)
    {
      biginteger modulus ;
      modulus.setValue(va.modulus[0].getValueTemp());
      for(unsigned int i = 1 ; i < va.modulus.size() ; ++i)
	{
	  if(modulus != va.modulus[i]) // if one is different: no modulus
	    return(bigintegerR::create_SEXP(result));
	}
      result.modulus.push_back(modulus);
    }



  return(bigintegerR::create_SEXP(result));

}

 
// Return min
SEXP biginteger_min(SEXP a, SEXP narm)
{
  bigvec result;

  bigvec va = bigintegerR::create_bignum(a);

  if( ! va.size())
    return bigintegerR::create_SEXP(result);

  unsigned int minimum = 0;

  PROTECT (a = AS_INTEGER(a));
  int na_remove = INTEGER(a)[0];
  UNPROTECT(1);


  for(unsigned int i = 1 ; i < va.size(); ++i)
    {
      if(va.value[i].isNA() && 
	 ( ! na_remove) )
	return bigintegerR::create_SEXP(result);
      else
	if(!(va.value[i] >  va.value[minimum] ))
	  minimum = i; 
    }
  
  result.push_back(va.value[minimum]);

  // now the modulus !
  if(va.modulus.size() == 1)
    result.modulus.push_back(va.modulus[0]);

  if(va.modulus.size()>1)
    {
      biginteger modulus ;
      modulus.setValue(va.modulus[0].getValueTemp());
      for(unsigned int i = 1 ; i < va.modulus.size() ; ++i)
	{
	  if(modulus != va.modulus[i]) // if one is different: no modulus
	    return bigintegerR::create_SEXP(result);
	}
      result.modulus.push_back(modulus);
    }

  return bigintegerR::create_SEXP(result);

}


SEXP biginteger_cumsum(SEXP a)
{

  bigvec result;

  bigvec va = bigintegerR::create_bignum(a);

  result.value.resize(va.value.size());


  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
 
  bool hasmodulus = true;

  // first the modulus !
  if(va.modulus.size()>1)
    {
      biginteger modulus ;
      modulus.setValue(va.modulus[0].getValueTemp());
      
      for(unsigned int i = 1 ; i < va.modulus.size() ; ++i)
	{
	  if(modulus != va.modulus[i]) // if one is different: no modulus
	    hasmodulus = false;
	}
      if(hasmodulus)
	result.modulus.push_back(modulus);
      
    }
  else
    hasmodulus = false;
  
  if(va.modulus.size() == 1)
    {
      result.modulus.push_back(va.modulus[0]);
      hasmodulus = true;
   }




  for(unsigned int i = 0 ; i < va.size(); ++i)
    {
      {
	if(va.value[i].isNA() )
	  {	
	    break; // all last values are NA.
	  }
      
	mpz_add(val,val,va.value[i].getValueTemp());
	
	if(hasmodulus)
	  mpz_mod(val,val,va.modulus[0].getValueTemp() );
	
	result.value[i].setValue(val);
      }
    }

  return(bigintegerR::create_SEXP(result));

}


SEXP biginteger_prod(SEXP a)
{

  bigvec result;

  bigvec va = bigintegerR::create_bignum(a);

  result.value.resize(1);


  mpz_t val;
  mpz_init(val);
  mpz_set_ui(val,1);
  mpz_t_sentry val_s(val);
 
  bool hasmodulus = true;

  // first the modulus !
  if(va.modulus.size()>1)
    {
      biginteger modulus ;
      modulus.setValue(va.modulus[0].getValueTemp());
      
      for(unsigned int i = 1 ; i < va.modulus.size() ; ++i)
	{
	  if(modulus != va.modulus[i]) // if one is different: no modulus
	    hasmodulus = false;
	}
      if(hasmodulus)
	result.modulus.push_back(modulus);
      
    }
  else
    hasmodulus = false;
  
  if(va.modulus.size() == 1)
    {
      result.modulus.push_back(va.modulus[0]);
      hasmodulus = true;
    }




  for(unsigned int i = 0 ; i < va.size(); ++i)
    {
      {
	if(va.value[i].isNA() )
	  {	
	    return (bigintegerR::create_SEXP(result));
	  }
      
	mpz_mul(val,val,va.value[i].getValueTemp());
	
	if(hasmodulus)
	  mpz_mod(val,val,va.modulus[0].getValueTemp() );
	
      }
    }

  result.value[0].setValue(val);

  return(bigintegerR::create_SEXP(result));

}



