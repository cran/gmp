/************************************************************/
/*! \file bigintegerR.cc
 *  \brief C function to interface R and libgmp with big integer values
 *
 *  \version 1
 *
 *  \date Created: 27/10/04   
 *  \date Last modified: Time-stamp: <2005-02-27 09:49:01 antoine>
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

#include "biginteger.h"


#include <stdio.h>

#include <vector>
#include <algorithm>
using namespace std;
#include "bigintegerR.h"

/* Globals variables */

gmp_randstate_t seed_state;
int seed_init=0;


extern "C"
{
    /**
     * \brief Addition of a and b
     */
    SEXP biginteger_add(SEXP a, SEXP b);

    /**
     * \brief Subtraction of a and b
     */
    SEXP biginteger_sub(SEXP a, SEXP b);

    /**
     * \brief Multiplication of a and b
     */
    SEXP biginteger_mul(SEXP a, SEXP b);

    /**
     * \brief Quotient of a / b
     */
    SEXP biginteger_div(SEXP a, SEXP b);

    /**
     * \brief Modulus of a % b
     */
    SEXP biginteger_mod(SEXP a, SEXP b);

    /**
     * \brief Power of base a to exponent b
     */
    SEXP biginteger_pow(SEXP a, SEXP b);

    /**
     * \brief Inverse from a mod b
     */
    SEXP biginteger_inv(SEXP a, SEXP b);

    /**
     * \brief Greatest common divisor of a and b
     */
    SEXP biginteger_gcd(SEXP a, SEXP b);

    /**
     * \brief Least common multiply of a and b
     */
    SEXP biginteger_lcm(SEXP a, SEXP b);
    
    /**
     * \brief Sets the intern modulus of a to b
     */
    SEXP biginteger_setmod(SEXP a, SEXP b);

   
    /**
     * \brief Return from vector a all elements specified in vector b
     */
    SEXP biginteger_get_at(SEXP a, SEXP b);
    
    /**
     * \brief Return a vector with the values from src specified by 
     * idx to sequentiell values from "value".
     */
    SEXP biginteger_set_at(SEXP src, SEXP idx, SEXP value);
    
    /**
     * \brief Convert from an long value or a string into biginteger.
     *
     * If you want a modulus-set-bigint, just use 
     * as.biginteger(value, modulus)
     */
    SEXP biginteger_as(SEXP a, SEXP mod);

    /**
     * \brief Convert from a biginteger vector to a character string vector.
     */
    SEXP biginteger_as_character(SEXP a,SEXP b);

    /**
     * \brief Convert from a biginteger vector to a real vector.
     */
    SEXP biginteger_as_numeric(SEXP a);

    /**
     * \brief Return the length of the vector
     */
    SEXP biginteger_length(SEXP a);

    /**
     * \brief Returns a resized vector cut at end or filled with NA.
     */
    SEXP biginteger_setlength(SEXP vec, SEXP value);

    /**
     * \brief Return whether the parameter is NA
     */
    SEXP biginteger_is_na(SEXP a);


    /**
     * \brief Return sign of a
     */
    SEXP biginteger_sgn(SEXP a) ;

    /**
     * \brief Return whether a < b
     */
    SEXP biginteger_lt(SEXP a, SEXP b);
    
    /**
     * \brief Return whether a > b
     */
    SEXP biginteger_gt(SEXP a, SEXP b);
    
    /**
     * \brief Return whether a <= b
     */
    SEXP biginteger_lte(SEXP a, SEXP b);
    
    /**
     * \brief Return whether a >= b
     */
    SEXP biginteger_gte(SEXP a, SEXP b);
    
    /**
     * \brief Return whether a == b
     */
    SEXP biginteger_eq(SEXP a, SEXP b);
    
    /**
     * \brief Return whether a != b
     */
    SEXP biginteger_neq(SEXP a, SEXP b);

    /**
     * \brief For function c()
     */
    SEXP biginteger_c(SEXP args) ;

    /**
     * \brief Create vector as n times x
     */
    SEXP biginteger_rep(SEXP x, SEXP times) ;

    /**
     * \brief Return if a is prime with a proba test
     */
    SEXP biginteger_is_prime(SEXP a, SEXP reps) ;  

    /**
     * \brief Return next prime with a proba test
     */
    SEXP biginteger_nextprime(SEXP a) ;


    /**
     * \brief Return absolute value of a
     */
    SEXP biginteger_abs(SEXP a) ;

    /**
     * \brief Return bezoult coefficients
     */
    SEXP biginteger_gcdex(SEXP a, SEXP b) ;

    /** 
     * \brief Random number generation
     */ 
    SEXP biginteger_rand_u (SEXP nb ,SEXP length,SEXP newseed, SEXP ok);

    /**
     * \brief Give size of integer 
     */
    SEXP biginteger_sizeinbase (SEXP x, SEXP exp);


    /**
     * \brief Compute Fibonacci number
     */
    SEXP fibnum(SEXP n) ;

    /**
     * \brief Compute Fibonacci number
     */

    SEXP fibnum2(SEXP n) ;
    /**
     * \brief Compute lucas number
     */
    SEXP lucnum(SEXP n) ;

    /**
     * \brief Compute lucas number
     */
      SEXP lucnum2(SEXP n) ;

}


namespace bigintegerR
{
  /** \brief create a vector of bigmods, all without a modulus.
   */
    vector<bigmod> create_vector(SEXP param) {
	switch (TYPEOF(param)) {
	case RAWSXP:
	    {
		// deserialise the vector. first int is the size.
		vector<bigmod> v;
		char* raw = (char*)RAW(param);
		int pos = sizeof(int); // position in raw[]. Starting after header.
		for (int i = 0; i < ((int*)raw)[0]; ++i) {
		    v.push_back(bigmod(biginteger((void*)&raw[pos])));
		    pos += v.back().value.raw_size(); // increment number of bytes read.
		}
		return v;
	    }
	case REALSXP:
	    {
		double* d = REAL(param);
		vector<bigmod> v(d,d+LENGTH(param));
		for (int j = 0; j < v.size(); ++j)
		    if (d[j] == NA_REAL)
			v[j].value.setValue();
		return v;
	    }
	case INTSXP:
	case LGLSXP:
	    {
		int* i = INTEGER(param);
		vector<bigmod> v(i,i+LENGTH(param));
		for (int j = 0; j < v.size(); ++j)
		    if (i[j] == NA_INTEGER)
			v[j].value.setValue();
		return v;
	    }
	case STRSXP:
	    {
		vector<bigmod> v;
		v.reserve(LENGTH(param));
		for (int i = 0; i < LENGTH(param); ++i) {
		    if (STRING_ELT(param,i) == NA_STRING)
			v.push_back(biginteger());
		    else
			v.push_back(biginteger(std::string(CHAR(STRING_ELT(param,i)))));
		}
		return v;
	    }
	}
    }

    vector<bigmod> create_bignum(SEXP param) {
	SEXP modName;
	PROTECT(modName = Rf_allocVector(STRSXP,1));
	SET_STRING_ELT(modName, 0, Rf_mkChar("mod"));
	SEXP modAttr = Rf_getAttrib(param, modName);
	UNPROTECT(1);
    	if (TYPEOF(modAttr) != NILSXP) {
	    vector<bigmod> v = bigintegerR::create_vector(param);
	    vector<bigmod> attrib = bigintegerR::create_vector(modAttr);
	    if (!attrib.empty()) // sanity check
		for (int i = 0; i < v.size(); ++i)
		    v[i] = set_modulus(v[i], attrib[i%attrib.size()]);
	    return v;
	} else
	    return bigintegerR::create_vector(param);
    }

    vector<int> create_int(SEXP param) {
	switch (TYPEOF(param)) {
	case REALSXP:
	    {
		double* d = REAL(param);
		// copy vector manually to avoid stupid conversation warning in STL-code :-/
		vector<int> v;
		v.reserve(LENGTH(param));
		for (int i = 0; i < LENGTH(param); ++i)
		    v.push_back((int)d[i]);
		return v;
		//return vector<int>(d, d+LENGTH(param));
	    }
	case INTSXP:
	case LGLSXP:
	    {
		int* i = INTEGER(param);
		return vector<int>(i, i+LENGTH(param));
	    }
	default:
	    return vector<int>();
	}
    }


    SEXP create_SEXP(const vector<bigmod>& v) {
	SEXP ans;
	int size = sizeof(int); // starting with vector-size-header
	for (int i = 0; i < v.size(); ++i)
	    size += v[i].value.raw_size(); // adding each bigint's needed size
	PROTECT(ans = Rf_allocVector(RAWSXP, size));
	char* r = (char*)RAW(ans);
	((int*)r)[0] = v.size(); // first int is vector-size-header
	int pos = sizeof(int); // current position in r[] (starting after vector-size-header)
	for (int i = 0; i < v.size(); ++i)
	    pos += v[i].value.as_raw(&r[pos]);
	
	// set the class attribute to "bigz"
	SEXP className;
	PROTECT(className = Rf_allocVector(STRSXP, 1));
	SET_STRING_ELT(className, 0, Rf_mkChar("bigz"));
	Rf_setAttrib(ans, R_ClassSymbol, className);

	// set the mod attribute
	vector<bigmod> attrib;
	attrib.reserve(v.size());
	bool allNA = true;
	for (vector<bigmod>::const_iterator it = v.begin(); it != v.end(); ++it) {
	    if (!it->modulus.isNA())
		allNA = false;
	    attrib.push_back(it->modulus);
	}

	if (!allNA) {
	    SEXP modAttr = create_SEXP(attrib); 
	    SEXP modName;
	    PROTECT(modName = Rf_allocVector(STRSXP,1));
	    SET_STRING_ELT(modName, 0, Rf_mkChar("mod"));
	    Rf_setAttrib(ans, modName, modAttr);
	    UNPROTECT(1);
	}
	
	UNPROTECT(2);
	return ans;
    }

    /// A pointer to a binary operator for bigintegers
    typedef bigmod (*biginteger_binary_fn)(const bigmod&, const bigmod&);
    /**
     * \brief Main function of doing a binary operation on bigintegers.
     * It calls a function argument for doing the correct thing.
     * This could also be written as a class functor (template)
     * to save one function call, but then code bloat will happen.
     */
    SEXP biginteger_binary_operation(SEXP a, SEXP b, biginteger_binary_fn f)
    {
	vector<bigmod> va = bigintegerR::create_bignum(a);
	vector<bigmod> vb = bigintegerR::create_bignum(b);
	if (va.empty() || vb.empty())
	    Rf_error("argument must not be an empty list.");
	vector<bigmod> result;
	int size = max(va.size(), vb.size());
	result.reserve(size);
	for (int i = 0; i < size; ++i)
	    result.push_back(f(va[i%va.size()], vb[i%vb.size()]));
	return bigintegerR::create_SEXP(result);
    }

    typedef bool (*biginteger_logical_binary_fn)(const bigmod&, const bigmod&);
    SEXP biginteger_logical_binary_operation(SEXP a, SEXP b, biginteger_logical_binary_fn f)
    {
	vector<bigmod> va = bigintegerR::create_bignum(a);
	vector<bigmod> vb = bigintegerR::create_bignum(b);
	if (va.empty() || vb.empty())
	    Rf_error("argument must not be an empty list.");
	int size = max(va.size(), vb.size());
	SEXP ans;
	PROTECT(ans = Rf_allocVector(LGLSXP, size));
	for (int i = 0; i < size; ++i) {
	    bigmod am = va[i % va.size()];
	    bigmod bm = vb[i % vb.size()];
	    if (am.value.isNA() || bm.value.isNA())
		LOGICAL(ans)[i] = NA_LOGICAL;
	    else
		LOGICAL(ans)[i] = f(va[i%va.size()], vb[i%vb.size()]) ? 1 : 0;
	}
	UNPROTECT(1);
	return ans;
    }

    bool lt(const bigmod& lhs, const bigmod& rhs)
    {return mpz_cmp(lhs.value.getValueTemp(), rhs.value.getValueTemp()) < 0;}
    bool gt(const bigmod& lhs, const bigmod& rhs) 
    {return mpz_cmp(lhs.value.getValueTemp(), rhs.value.getValueTemp()) > 0;}
    bool lte(const bigmod& lhs, const bigmod& rhs)
    {return mpz_cmp(lhs.value.getValueTemp(), rhs.value.getValueTemp()) <= 0;}
    bool gte(const bigmod& lhs, const bigmod& rhs)
    {return mpz_cmp(lhs.value.getValueTemp(), rhs.value.getValueTemp()) >= 0;}
    bool eq(const bigmod& lhs, const bigmod& rhs)
    {return mpz_cmp(lhs.value.getValueTemp(), rhs.value.getValueTemp()) == 0;}
    bool neq(const bigmod& lhs, const bigmod& rhs)
    {return mpz_cmp(lhs.value.getValueTemp(), rhs.value.getValueTemp()) != 0;}

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
    vector<bigmod> v = bigintegerR::create_bignum(a);
    SEXP ans;
    int base;
    base = *INTEGER(b);
    PROTECT(ans = Rf_allocVector(STRSXP, v.size()));
    for (int i = 0; i < v.size(); ++i)
	SET_STRING_ELT(ans, i, Rf_mkChar(v[i].str(base).c_str()));
    
    UNPROTECT(1);
    return ans;
}

SEXP biginteger_as_numeric(SEXP a)
{
    vector<bigmod> v = bigintegerR::create_bignum(a);
    SEXP ans;
    PROTECT(ans = Rf_allocVector(REALSXP,v.size()));
    for (int i = 0; i < v.size(); ++i)
	REAL(ans)[i] = v[i].value.as_double();
    UNPROTECT(1);
    return ans;
}

SEXP biginteger_get_at(SEXP a, SEXP b)
{
    vector<bigmod> va = bigintegerR::create_bignum(a);
    vector<int> vb = bigintegerR::create_int(b);
    vector<bigmod> result;
    if (TYPEOF(b) == LGLSXP) {
	for (int i = 0; i < va.size(); ++i)
	    if (vb[i%vb.size()])
		result.push_back(va[i]);
    } else {
	vb.erase(remove(vb.begin(), vb.end(), 0), vb.end()); // remove all zeroes
	if (vb.empty())
	    return bigintegerR::create_SEXP(vector<bigmod>());
	if (vb[0] < 0) {
	    for (vector<int>::iterator it = vb.begin(); it != vb.end(); ++it)
		if (*it > 0)
		    Rf_error("only 0's may mix with negative subscripts");
		else if (-(*it)-1 >= va.size())
		    Rf_error("subscript out of bounds");
	    // TODO: This is optimized for large va.size and small vb.size.
	    // Maybe add a condition to use a different approach for large vb's
	    result.reserve(va.size()-vb.size());
	    for (int i = 0; i < va.size(); ++i)
		if (find(vb.begin(), vb.end(), -i-1) == vb.end())
		    result.push_back(va[i]);
	} else {
	    result.reserve(vb.size());
	    for (vector<int>::iterator it = vb.begin(); it != vb.end(); ++it) {
		if (*it < 0)
		    Rf_error("only 0's may mix with negative subscripts");
		if (*it <= va.size())
		    result.push_back(va[(*it)-1]);
		else
		    result.push_back(bigmod()); // NA for out of range's
	    }
	}
    }
    return bigintegerR::create_SEXP(result);
}

SEXP biginteger_set_at(SEXP src, SEXP idx, SEXP value)
{
    vector<bigmod> result = bigintegerR::create_bignum(src);
    vector<bigmod> vvalue = bigintegerR::create_bignum(value);
    vector<int> vidx = bigintegerR::create_int(idx);
    if (TYPEOF(idx) == LGLSXP) {
	int pos = 0;
	for (int i = 0; i < result.size(); ++i)
	    if (vidx[i%vidx.size()])
		result[i] = vvalue[pos++%vvalue.size()];
    } else {
	vidx.erase(remove(vidx.begin(), vidx.end(), 0), vidx.end()); // remove all zeroes
	if (vidx.empty())
	    return bigintegerR::create_SEXP(result);
	if (vidx[0] < 0) {
	    for (vector<int>::iterator it = vidx.begin(); it != vidx.end(); ++it)
		if (*it > 0)
		    Rf_error("only 0's may mix with negative subscripts");
		else if (-(*it)-1 >= result.size())
		    Rf_error("subscript out of bounds");
	    int pos = 0;
	    for (int i = 0; i < result.size(); ++i)
		if (find(vidx.begin(), vidx.end(), -i-1) == vidx.end())
		    result[i] = vvalue[pos++%vvalue.size()];
	} else {
	    // finding maximum to resize vector if needed
	    int maximum = INT_MIN;
	    for (vector<int>::iterator it = vidx.begin(); it != vidx.end(); ++it)
		maximum = max(maximum, *it);
	    if (maximum > result.size())
		result.resize(maximum);
	    int pos = 0;
	    for (vector<int>::iterator it = vidx.begin(); it != vidx.end(); ++it) {
		if (*it < 0)
		    Rf_error("only 0's may mix with negative subscripts");
		result[(*it)-1] = vvalue[pos++%vvalue.size()];
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
	    else if (len == NA_REAL)
		Rf_error("vector size cannot be NA");
	    break;
	case STRSXP:
	    // dunno why R spits out this strange error on "length(foo) <- -1"
	    // but I always follow the holy standard ;-)
	    Rf_error("negative length vectors are not allowed");
	default:
	    Rf_error("invalid second argument");
    }
    vector<bigmod> v =bigintegerR::create_bignum(vec);
    v.resize(len);
    return bigintegerR::create_SEXP(v);
}

SEXP biginteger_is_na(SEXP a) 
{
    vector<bigmod> v = bigintegerR::create_bignum(a);
    SEXP ans;
    PROTECT(ans = Rf_allocVector(LGLSXP, v.size()));
    for (int i = 0; i < v.size(); ++i)
	LOGICAL(ans)[i] = v[i].value.isNA();
    UNPROTECT(1);
    return ans;
}


SEXP biginteger_sgn(SEXP a) 
{
    vector<bigmod> v = bigintegerR::create_bignum(a);
    SEXP ans;
    PROTECT(ans = Rf_allocVector(INTSXP, v.size()));
    for (int i = 0; i < v.size(); ++i)
	INTEGER(ans)[i] = mpz_sgn(v[i].value.getValueTemp());
    UNPROTECT(1);
    return ans;
}

SEXP biginteger_c(SEXP args) 
{
    int i=0,j=0,size=0; 
    vector<bigmod> result;
    vector<bigmod> v;

    /* Count number of all elements */
    for(i = 0; i< LENGTH(args); i++)
	size += LENGTH(VECTOR_ELT(args,i));

    result.reserve(size);

    for(i =0; i<LENGTH(args);i++)
      {
	v = bigintegerR::create_bignum(VECTOR_ELT(args,i));
	for(j=0; j< v.size() ; j++)
	    result.push_back(v[j]);
	v.clear();
      }

    return bigintegerR::create_SEXP(result);
}


SEXP biginteger_rep(SEXP x, SEXP times) 
{
  vector<bigmod> v = bigintegerR::create_bignum(x);
  vector<bigmod> result;
  int i,j,size,rep;
  
  PROTECT(times= AS_INTEGER(times));
  rep = INTEGER(times)[0];
  UNPROTECT(1);

  result.reserve(v.size()*rep);
  for(i = 0 ; i< rep ; i++)
    for(j = 0 ; j < v.size() ; j++)
      result.push_back(v[j]);
  
  return bigintegerR::create_SEXP(result);
}



SEXP biginteger_is_prime(SEXP a, SEXP reps) 
{
    vector<bigmod> v = bigintegerR::create_bignum(a);
    vector<int> vb = bigintegerR::create_int(reps);
    SEXP ans;
    PROTECT(ans = Rf_allocVector(INTSXP, v.size()));
    if(v.size() == vb.size())
      for (int i = 0; i < v.size(); ++i)
	INTEGER(ans)[i] = v[i].value.isprime(vb[i]);
    else
      for (int i = 0; i < v.size(); ++i)
	INTEGER(ans)[i] = v[i].value.isprime(vb[0]);
    UNPROTECT(1);
    return ans;
}
 
/** brief Use this to clear mpz_t structs at end-of-function automatically
 */
struct mpz_t_sentry {
  mpz_t& value;
  mpz_t_sentry(mpz_t& v): value(v) {}
  ~mpz_t_sentry() {mpz_clear(value);}
};
 

SEXP biginteger_nextprime(SEXP a) 
{
    vector<bigmod> v =bigintegerR::create_bignum(a);
    vector<bigmod> result;
    
    result.reserve(v.size());
    
    mpz_t val;
    mpz_init(val);
    mpz_t_sentry val_s(val);

    for (int i = 0; i < v.size(); ++i)
      {
	mpz_nextprime(val,v[i].value.getValueTemp());
	result.push_back(bigmod());
	result[i].value.setValue(val);
      }
    
    return bigintegerR::create_SEXP(result);
    
}

SEXP biginteger_abs(SEXP a) 
{
    vector<bigmod> v =bigintegerR::create_bignum(a);
    vector<bigmod> result;
    
    result.reserve(v.size());
    
    mpz_t val;
    mpz_init(val);
    mpz_t_sentry val_s(val);

    for (int i = 0; i < v.size(); ++i)
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
    vector<bigmod> va = bigintegerR::create_bignum(a);
    vector<bigmod> vb = bigintegerR::create_bignum(b);
    vector<bigmod> result;
    int i;

    if(va.size() != vb.size())
      return bigintegerR::create_SEXP(vector<bigmod>());
    

    result.reserve(3*va.size());

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
	result.push_back(bigmod()); // Hem... not very elegant !
	result.push_back(bigmod());
	result.push_back(bigmod());
	result[i*3].value.setValue(g);
	result[i*3+1].value.setValue(s);
	result[i*3+2].value.setValue(t);
	
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
  SEXP ans;

  mpz_t    bz;
  char * mot;
  int i,l,flag,len,size;
  vector<bigmod> result;
   
  extern int seed_init;
  extern gmp_randstate_t seed_state;


  /* store input data into appropriate mode */
  vector<bigmod> va = bigintegerR::create_bignum(newseed);
  PROTECT (ok = AS_INTEGER(ok));
  PROTECT (length = AS_INTEGER(length));
  PROTECT (nb = AS_INTEGER(nb));
  flag = INTEGER(ok)[0];
  len = INTEGER(length)[0];
  size = INTEGER(nb)[0];
  UNPROTECT(3);

  result.reserve(size);
  
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
      result.push_back(bigmod());
      result[i].value.setValue(bz);
    }
    return bigintegerR::create_SEXP(result);
}


/** @brief biginteger_sizeinbase return 
 *  @param x BigInteger
 *  @param base BigInteger
 */
SEXP biginteger_sizeinbase(SEXP x, SEXP base) 
{
    vector<bigmod> vx = bigintegerR::create_bignum(x);

    int i,basesize;
    double reslog;
    size_t taille;

    SEXP ans;
    
    PROTECT (base = AS_INTEGER(base));
    basesize=INTEGER(base)[0];
    UNPROTECT(1);

    PROTECT(ans = Rf_allocVector(INTSXP,vx.size()));

    for(i=0; i<vx.size();i++)
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
  vector<bigmod> result;
  int num;

  PROTECT (n = AS_INTEGER(n));
  num = INTEGER(n)[0];
  UNPROTECT(1);

  result.reserve(1);

  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
  
  mpz_fib_ui(val,num);
  result.push_back(bigmod());
  result[0].value.setValue(val);


  return bigintegerR::create_SEXP(result);

}

/** @brief fibnum2 return nth and n-1th Fibonacci number 
 *  @param n integer
 */
SEXP fibnum2(SEXP n) 
{
  vector<bigmod> result;
  int num;

  PROTECT (n = AS_INTEGER(n));
  num = INTEGER(n)[0];
  UNPROTECT(1);

  result.reserve(1);

  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
  mpz_t val2;
  mpz_init(val2);
  mpz_t_sentry val_s2(val2);
  
  mpz_fib2_ui(val,val2,num);
  result.push_back(bigmod());
  result[0].value.setValue(val2);
  result.push_back(bigmod());
  result[1].value.setValue(val);


  return bigintegerR::create_SEXP(result);

}

/** @brief lucnum return nth lucas number 
 *  @param n integer
 */
SEXP lucnum(SEXP n) 
{
  vector<bigmod> result;
  int num;

  PROTECT (n = AS_INTEGER(n));
  num = INTEGER(n)[0];
  UNPROTECT(1);

  result.reserve(1);

  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
  
  mpz_lucnum_ui(val,num);
  result.push_back(bigmod());
  result[0].value.setValue(val);


  return bigintegerR::create_SEXP(result);

}

/** @brief lucnum2 return nth and n-1th lucas number 
 *  @param n integer
 */
SEXP lucnum2(SEXP n) 
{
  vector<bigmod> result;
  int num;

  PROTECT (n = AS_INTEGER(n));
  num = INTEGER(n)[0];
  UNPROTECT(1);

  result.reserve(1);

  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
  mpz_t val2;
  mpz_init(val2);
  mpz_t_sentry val_s2(val2);
  
  mpz_lucnum2_ui(val,val2,num);
  result.push_back(bigmod());
  result[0].value.setValue(val2);
  result.push_back(bigmod());
  result[1].value.setValue(val);


  return bigintegerR::create_SEXP(result);

}


