/************************************************************/
/*! \file bigrationalR.cc
 *  \brief C function to interface R and libgmp with big rational values
 *
 *  \version 1
 *
 *  \date Created: 12/12/04   
 *  \date Last modified: Time-stamp: <2005-02-27 09:55:12 antoine>
 *
 *  \author Antoine Lucas (adapted from biginteger class made by
 *                         Immanuel Scholz)
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

#include "bigrational.h"


#include <stdio.h>

#include <vector>
#include <algorithm>
using namespace std;
#include "bigintegerR.h"


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

    /**
     * \brief Create vector as n times x
     */
    SEXP bigrational_rep(SEXP x, SEXP times) ;

}


namespace bigrationalR
{
  /** \brief create a vector of bigrationals, all without a denominator.
   */
    vector<bigrational> create_vector(SEXP param) {
	switch (TYPEOF(param)) {
	case RAWSXP:
	    {
		// deserialise the vector. first int is the size.
		vector<bigrational> v;
		char* raw = (char*)RAW(param);
		int pos = sizeof(int); // position in raw[]. Starting after header.
		for (int i = 0; i < ((int*)raw)[0]; ++i) {
		    v.push_back(bigrational(bigrational((void*)&raw[pos])));
		    pos += v.back().numerator.raw_size(); // increment number of bytes read.
		}
		return v;
	    }
	case REALSXP:
	    {
		double* d = REAL(param);
		vector<bigrational> v(d,d+LENGTH(param));
		for (int j = 0; j < v.size(); ++j)
		    if (d[j] == NA_REAL)
			v[j].numerator.setValue();
		    else
		      v[j].setQValue(d[j]);
		return v;
	    }
	case INTSXP:
	case LGLSXP:
	    {
		int* i = INTEGER(param);
		vector<bigrational> v(i,i+LENGTH(param));
		for (int j = 0; j < v.size(); ++j)
		    if (i[j] == NA_INTEGER)
			v[j].numerator.setValue();
		return v;
	    }
	case STRSXP:
	    {
		vector<bigrational> v;
		v.reserve(LENGTH(param));
		for (int i = 0; i < LENGTH(param); ++i) {
		    if (STRING_ELT(param,i) == NA_STRING)
			v.push_back(bigrational());
		    else
			v.push_back(bigrational(std::string(CHAR(STRING_ELT(param,i)))));
		}
		return v;
	    }
	}
    }

    vector<bigrational> create_bignum(SEXP param) {
	SEXP denName;
	PROTECT(denName = Rf_allocVector(STRSXP,1));
	SET_STRING_ELT(denName, 0, Rf_mkChar("denominator"));
	SEXP denAttr = Rf_getAttrib(param, denName);
	UNPROTECT(1);
    	if (TYPEOF(denAttr) != NILSXP) {
	    vector<bigrational> v = bigrationalR::create_vector(param);
	    vector<bigrational> attrib = bigrationalR::create_vector(denAttr);
	    if (!attrib.empty()) // sanity check
		for (int i = 0; i < v.size(); ++i)
		  {
		    /*v[i] = set_denominator(v[i], attrib[i%attrib.size()]);*/
		    if(v[i].denominator.isNA() && (attrib[i%attrib.size()].numerator.sgn() != 0) )
		      v[i].denominator.setValue ((__mpz_struct*)attrib[i%attrib.size()].numerator.getValueTemp()) ; 
		    /*printf("q:%f\n",v[i].as_double());*/
		  }
	    return v;
	} else
	    return bigrationalR::create_vector(param);
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


    SEXP create_SEXP(const vector<bigrational>& v) {
	SEXP ans;
	int size = sizeof(int); // starting with vector-size-header
	for (int i = 0; i < v.size(); ++i)
	    size += v[i].numerator.raw_size(); // adding each bigint's needed size
	PROTECT(ans = Rf_allocVector(RAWSXP, size));
	char* r = (char*)RAW(ans);
	((int*)r)[0] = v.size(); // first int is vector-size-header
	int pos = sizeof(int); // current position in r[] (starting after vector-size-header)
	for (int i = 0; i < v.size(); ++i)
	    pos += v[i].numerator.as_raw(&r[pos]);
	
	// set the class attribute to "bigrational"
	SEXP className;
	PROTECT(className = Rf_allocVector(STRSXP, 1));
	SET_STRING_ELT(className, 0, Rf_mkChar("bigq"));
	Rf_setAttrib(ans, R_ClassSymbol, className);

	// set the denominator attribute
	vector<bigrational> attrib;
	attrib.reserve(v.size());
	bool allNA = true;
	for (vector<bigrational>::const_iterator it = v.begin(); it != v.end(); ++it) {
	    if (!it->denominator.isNA())
		allNA = false;
	    attrib.push_back(it->denominator);
	}

	if (!allNA) {
	    SEXP denAttr = create_SEXP(attrib); 
	    SEXP denName;
	    PROTECT(denName = Rf_allocVector(STRSXP,1));
	    SET_STRING_ELT(denName, 0, Rf_mkChar("denominator"));
	    Rf_setAttrib(ans, denName, denAttr);
	    UNPROTECT(1);
	}
	
	UNPROTECT(2);
	return ans;
    }

    /// A pointer to a binary operator for bigrationals
    typedef bigrational (*bigrational_binary_fn)(const bigrational&, const bigrational&);
    /**
     * \brief Main function of doing a binary operation on bigrationals.
     * It calls a function argument for doing the correct thing.
     * This could also be written as a class functor (template)
     * to save one function call, but then code bloat will happen.
     */
    SEXP bigrational_binary_operation(SEXP a, SEXP b, bigrational_binary_fn f)
    {
	vector<bigrational> va = bigrationalR::create_bignum(a);
	vector<bigrational> vb = bigrationalR::create_bignum(b);
	if (va.empty() || vb.empty())
	    Rf_error("argument must not be an empty list.");
	vector<bigrational> result;
	int size = max(va.size(), vb.size());
	result.reserve(size);
	for (int i = 0; i < size; ++i)
	    result.push_back(f(va[i%va.size()], vb[i%vb.size()]));
	return bigrationalR::create_SEXP(result);
    }

    typedef bool (*bigrational_logical_binary_fn)(const bigrational&, const bigrational&);
    SEXP bigrational_logical_binary_operation(SEXP a, SEXP b, bigrational_logical_binary_fn f)
    {
	vector<bigrational> va = bigrationalR::create_bignum(a);
	vector<bigrational> vb = bigrationalR::create_bignum(b);
	if (va.empty() || vb.empty())
	    Rf_error("argument must not be an empty list.");
	int size = max(va.size(), vb.size());
	SEXP ans;
	PROTECT(ans = Rf_allocVector(LGLSXP, size));
	for (int i = 0; i < size; ++i) {
	    bigrational am = va[i % va.size()];
	    bigrational bm = vb[i % vb.size()];
	    if (am.numerator.isNA() || bm.numerator.isNA())
		LOGICAL(ans)[i] = NA_LOGICAL;
	    else
		LOGICAL(ans)[i] = f(va[i%va.size()], vb[i%vb.size()]) ? 1 : 0;
	}
	UNPROTECT(1);
	return ans;
    }

    bool lt(const bigrational& lhs, const bigrational& rhs)
    {
      int res;
      mpq_t lhsq, rhsq;
      GET_MPQ (lhs,lhsq,lhsqt)
      GET_MPQ (rhs,rhsq,rhsqt)
      return mpq_cmp(lhsq, rhsq) < 0;
    }
    bool gt(const bigrational& lhs, const bigrational& rhs) 
    {
      int res;
      mpq_t lhsq, rhsq;
      GET_MPQ (lhs,lhsq,lhsqt)
      GET_MPQ (rhs,rhsq,rhsqt)
      return mpq_cmp(lhsq, rhsq) > 0;
    }
    bool lte(const bigrational& lhs, const bigrational& rhs)
    {
      int res;
      mpq_t lhsq, rhsq;
      GET_MPQ (lhs,lhsq,lhsqt)
      GET_MPQ (rhs,rhsq,rhsqt)
      return mpq_cmp(lhsq, rhsq) <= 0;
    }
    bool gte(const bigrational& lhs, const bigrational& rhs)
    {
      int res;
      mpq_t lhsq, rhsq;
      GET_MPQ (lhs,lhsq,lhsqt)
      GET_MPQ (rhs,rhsq,rhsqt)
      return mpq_cmp(lhsq, rhsq) >= 0;
    }
    bool eq(const bigrational& lhs, const bigrational& rhs)
    {
      int res;
      mpq_t lhsq, rhsq;
      GET_MPQ (lhs,lhsq,lhsqt)
      GET_MPQ (rhs,rhsq,rhsqt)
      return mpq_equal(lhsq, rhsq) != 0;
    }
    bool neq(const bigrational& lhs, const bigrational& rhs)
    {
      int res;
      mpq_t lhsq, rhsq;
      GET_MPQ (lhs,lhsq,lhsqt)
      GET_MPQ (rhs,rhsq,rhsqt)
      return mpq_equal(lhsq, rhsq) == 0;
    }


}

/* End of namespace bigrationalR*/



SEXP bigrational_add (SEXP a, SEXP b) {return bigrationalR::bigrational_binary_operation(a,b,operator+);}
SEXP bigrational_sub (SEXP a, SEXP b) {return bigrationalR::bigrational_binary_operation(a,b,operator-);}
SEXP bigrational_mul (SEXP a, SEXP b) {return bigrationalR::bigrational_binary_operation(a,b,operator*);}
SEXP bigrational_div (SEXP a, SEXP b) {return bigrationalR::bigrational_binary_operation(a,b,operator/);}
SEXP bigrational_as (SEXP n, SEXP d) {return bigrationalR::bigrational_binary_operation(n,d,set_denominator);}

SEXP bigrational_lt (SEXP a, SEXP b) {return bigrationalR::bigrational_logical_binary_operation(a,b,bigrationalR::lt);}
SEXP bigrational_gt (SEXP a, SEXP b) {return bigrationalR::bigrational_logical_binary_operation(a,b,bigrationalR::gt);}
SEXP bigrational_lte (SEXP a, SEXP b) {return bigrationalR::bigrational_logical_binary_operation(a,b,bigrationalR::lte);}
SEXP bigrational_gte (SEXP a, SEXP b) {return bigrationalR::bigrational_logical_binary_operation(a,b,bigrationalR::gte);}
SEXP bigrational_eq (SEXP a, SEXP b) {return bigrationalR::bigrational_logical_binary_operation(a,b,bigrationalR::eq);}
SEXP bigrational_neq (SEXP a, SEXP b) {return bigrationalR::bigrational_logical_binary_operation(a,b,bigrationalR::neq);}



SEXP bigrational_as_character(SEXP a,SEXP b)
{
    vector<bigrational> v = bigrationalR::create_bignum(a);
    SEXP ans;
    int base;
    base = *INTEGER(b);
    PROTECT(ans = Rf_allocVector(STRSXP, v.size()));
    for (int i = 0; i < v.size(); ++i)
	SET_STRING_ELT(ans, i, Rf_mkChar(v[i].str(base).c_str()));
    UNPROTECT(1);
    return ans;
}

SEXP bigrational_as_numeric(SEXP a)
{
    vector<bigrational> v = bigrationalR::create_bignum(a);
    SEXP ans;
    PROTECT(ans = Rf_allocVector(REALSXP,v.size()));
    for (int i = 0; i < v.size(); ++i)
      {
	mpq_t val;
	mpq_init(val);
	mpq_t_sentry val_s(val);
	mpq_set_num(val,v[i].numerator.getValueTemp());
	if(!v[i].denominator.isNA())
	  mpq_set_den(val,v[i].denominator.getValueTemp());

	REAL(ans)[i] = mpq_get_d(val);
      }
    UNPROTECT(1);
    return ans;
}

SEXP bigrational_get_at(SEXP a, SEXP b)
{
    vector<bigrational> va = bigrationalR::create_bignum(a);
    vector<int> vb = bigrationalR::create_int(b);
    vector<bigrational> result;
    if (TYPEOF(b) == LGLSXP) {
	for (int i = 0; i < va.size(); ++i)
	    if (vb[i%vb.size()])
		result.push_back(va[i]);
    } else {
	vb.erase(remove(vb.begin(), vb.end(), 0), vb.end()); // remove all zeroes
	if (vb.empty())
	    return bigrationalR::create_SEXP(vector<bigrational>());
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
		    result.push_back(bigrational()); // NA for out of range's
	    }
	}
    }
    return bigrationalR::create_SEXP(result);
}

SEXP bigrational_set_at(SEXP src, SEXP idx, SEXP value)
{
    vector<bigrational> result = bigrationalR::create_bignum(src);
    vector<bigrational> vvalue = bigrationalR::create_bignum(value);
    vector<int> vidx = bigrationalR::create_int(idx);
    if (TYPEOF(idx) == LGLSXP) {
	int pos = 0;
	for (int i = 0; i < result.size(); ++i)
	    if (vidx[i%vidx.size()])
		result[i] = vvalue[pos++%vvalue.size()];
    } else {
	vidx.erase(remove(vidx.begin(), vidx.end(), 0), vidx.end()); // remove all zeroes
	if (vidx.empty())
	    return bigrationalR::create_SEXP(result);
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
    return bigrationalR::create_SEXP(result);
}

SEXP bigrational_length(SEXP a)
{
    SEXP ans;
    PROTECT(ans = Rf_allocVector(INTSXP,1));
    INTEGER(ans)[0] = bigrationalR::create_bignum(a).size();
    UNPROTECT(1);
    return ans;
}

SEXP bigrational_den(SEXP a)
{
    vector<bigrational> v =bigrationalR::create_bignum(a);
    vector<bigmod> result;
    
    result.reserve(v.size());
    
    for (int i = 0; i < v.size(); ++i)
      {
	result.push_back(bigmod());
	result[i].value.setValue((__mpz_struct*)v[i].denominator.getValueTemp());
      }
    
    return bigintegerR::create_SEXP(result);
}



SEXP bigrational_num(SEXP a)
{
    vector<bigrational> v =bigrationalR::create_bignum(a);
    vector<bigmod> result;
    
    result.reserve(v.size());
    
    for (int i = 0; i < v.size(); ++i)
      {
	result.push_back(bigmod());
	result[i].value.setValue((__mpz_struct*)v[i].numerator.getValueTemp());
      }
    
    return bigintegerR::create_SEXP(result);
}

SEXP bigrational_setlength(SEXP vec, SEXP value)
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
    vector<bigrational> v =bigrationalR::create_bignum(vec);
    v.resize(len);
    return bigrationalR::create_SEXP(v);
}

SEXP bigrational_is_na(SEXP a) 
{
    vector<bigrational> v = bigrationalR::create_bignum(a);
    SEXP ans;
    PROTECT(ans = Rf_allocVector(LGLSXP, v.size()));
    for (int i = 0; i < v.size(); ++i)
	LOGICAL(ans)[i] = v[i].numerator.isNA();
    UNPROTECT(1);
    return ans;
}

SEXP bigrational_c(SEXP args) 
{
    int i=0,j=0,size=0; 
    vector<bigrational> result;
    vector<bigrational> v;

    /* Count number of all elements */
    for(i = 0; i< LENGTH(args); i++)
	size += LENGTH(VECTOR_ELT(args,i));

    result.reserve(size);

    for(i =0; i<LENGTH(args);i++)
      {
	v = bigrationalR::create_bignum(VECTOR_ELT(args,i));
	for(j=0; j< v.size() ; j++)
	    result.push_back(v[j]);
	v.clear();
      }

    return bigrationalR::create_SEXP(result);
}


SEXP bigrational_rep(SEXP x, SEXP times) 
{
  vector<bigrational> v = bigrationalR::create_bignum(x);
  vector<bigrational> result;
  int i,j,size,rep;
  
  PROTECT(times= AS_INTEGER(times));
  rep = INTEGER(times)[0];
  UNPROTECT(1);

  result.reserve(v.size()*rep);
  for(i = 0 ; i< rep ; i++)
    for(j = 0 ; j < v.size() ; j++)
      result.push_back(v[j]);
  
  return bigrationalR::create_SEXP(result);
}
