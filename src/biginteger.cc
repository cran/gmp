/*! \file biginteger.cc
 *  \brief C function for class biginteger & bigmod
 *
 *  \version 1
 *
 *  \date Created: 27/10/04   
 *  \date Last modified: Time-stamp: <2004-11-27 18:40:56 antoine>
 *
 *  \author Immanuel Scholz 
 *
 *  \note Licence: GPL
 */

#define USE_RINTERNALS
#define R_NO_REMAP   // avoid collisions with stl definitions

#include "biginteger.h"
#include <Rinternals.h>

#include <stdio.h>

using std::string;

biginteger::biginteger(void* raw)
{
    mpz_init(value);
    int* r = (int*)raw;
    if (r[0]) {
        mpz_import(value, r[0], 1, sizeof(int), 0, 0, &r[1]);
	na = false;
    } else
	setValue();
}

string biginteger::str() const
{
    if (isNA())
	return "NA";
    
    char* buf = new char[mpz_sizeinbase(value, 10)+2];
    mpz_get_str(buf, 10, value);
    string s = buf;
    delete [] buf;
    return s;
}

int biginteger::as_raw(void* raw) const
{
    int totals = raw_size();
    memset(raw, 0, totals);
    int* r = (int*)raw;
    r[0] = totals/sizeof(int) - 1;
    if (!isNA())
        mpz_export(&r[1], 0, 1, sizeof(int), 0, 0, value);
    return totals;
}

size_t biginteger::raw_size() const
{
    if (isNA())
	return sizeof(int);

    int numb = 8*sizeof(int);
    return sizeof(int) * (1 + (mpz_sizeinbase(value,2)+numb-1) / numb);
}

void biginteger::swap(biginteger& other)
{
    mpz_swap(value, other.value);
    bool n = na;
    na = other.na;
    other.na = n;
}


string bigmod::str() const
{
    if (value.isNA())
	return "NA";

    string s; // sstream seems to collide with libgmp :-(
    if (!modulus.isNA())
	s = "(";
    s += value.str();
    if (!modulus.isNA()) {
	s += " %% ";
	s += modulus.str();
	s += ")";
    }
    return s;
}

