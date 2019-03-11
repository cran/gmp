

#include "bigmod.h"

#define DEBUG_bigmod
#undef DEBUG_bigmod

#ifdef DEBUG_bigmod
#include <R_ext/Print.h>
#endif

// for GetOption etc:
#include <R.h>
#include <Rinternals.h>


// string representation of (.) wrt base 'b' :
std::string bigmod::str(int b) const
{
  if (value.isNA())
    return "NA";

  std::string s; // sstream seems to collide with libgmp :-(
  if (!modulus.isNA())
    s = "(";
  s += value.str(b);
  if (!modulus.isNA()) {
    s += " %% ";
    s += modulus.str(b);
    s += ")";
  }
  return s;
}

bigmod & bigmod::operator= (const bigmod& rhs)
{
  if(this != &rhs)
    {
      modulus.setValue( rhs.getModulus() );
      value.setValue(rhs.value );
    }
  return(*this);
}

bigmod bigmod::inv () const
{
  if(value.isNA() || modulus.isNA()) {
    return bigmod();
  }
    
  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
  if (mpz_invert(val, getValue().getValueTemp(), getModulus().getValueTemp()) == 0) {
    SEXP wOpt = Rf_GetOption1(Rf_install("gmp:warnNoInv"));
    if(wOpt != R_NilValue && Rf_asInteger(wOpt))
      warning(_("inv(x) returning NA as x has no inverse"));
    
    return bigmod(); // return NA; was
  }

  return bigmod(val, modulus);
}


bool operator!=(const bigmod& rhs, const bigmod& lhs)
{
  if(rhs.getValue() != lhs.getValue())
    return(true);
  return(rhs.getModulus() != lhs.getModulus());
}

bool operator==(const bigmod& rhs, const bigmod& lhs)
{
  if(rhs.getValue() != lhs.getValue())
    return(false);
  return(!(rhs.getModulus() != lhs.getModulus()));
}


DefaultBigMod operator+(const bigmod& lhs, const bigmod& rhs)
{
  return create_bigmod(lhs, rhs, mpz_add);
}

DefaultBigMod operator-(const bigmod& lhs, const bigmod& rhs)
{
  return create_bigmod(lhs, rhs, mpz_sub);
}

DefaultBigMod operator*(const bigmod& lhs, const bigmod& rhs)
{
  return create_bigmod(lhs, rhs, mpz_mul);
}

/* called via biginteger_binary_operation(.) from biginteger_div() in
 * ./bigintegerR.cc  from R's  .Call(biginteger_div,  a, b)
 *   ~~~~~~~~~~~~~~
 * itself called from  "/.bigz" = div.bigz()
 */
DefaultBigMod div_via_inv(const bigmod& a, const bigmod& b) {
    // compute  a/b  as  a * b^(-1)
    return operator*(a, pow(b, DefaultBigMod(-1)));
}


void integer_div(mpz_t result,const mpz_t a, const mpz_t b) {
  mpz_tdiv_q(result,a,b);
  //
  // si resulat < 0 et module != 0: on enleve 1.
  // i.e. 3 / -4 = -1
  if (mpz_sgn(a) * mpz_sgn(b) == -1) {
    mpz_t val;
    mpz_init(val);
    mpz_t_sentry val_s(val);
    mpz_mod(val, a, b);
    if (mpz_cmp_ui(val, 0) != 0) 
      {	
	mpz_sub_ui(result, result, 1);
      }
  }
}


/* called via biginteger_binary_operation(.) from R's
 * .Call(biginteger_divq, a, b) , itself called from '%/%.bigz' = divq.bigz()
 */
DefaultBigMod operator/(const bigmod& lhs, const bigmod& rhs) {
  return create_bigmod(lhs, rhs, integer_div, false);
}

DefaultBigMod operator%(const bigmod& lhs, const bigmod& rhs)
{
  if (lhs.getValue().isNA() || rhs.getValue().isNA())
    return DefaultBigMod();
  if (mpz_sgn(rhs.getValue().getValueTemp()) == 0) {
    warning(_("biginteger division by zero: returning NA"));
    return DefaultBigMod();
  }
  biginteger mod;
  if (!lhs.getModulus().isNA() || !rhs.getModulus().isNA())
    mod = rhs.getValue();

  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
  mpz_mod(val, lhs.getValue().getValueTemp(), rhs.getValue().getValueTemp());
  return DefaultBigMod(val, mod);
}


// Either 'base' has a modulus, or it has not *and* exp >= 0 :

DefaultBigMod pow(const bigmod& base, const bigmod& exp)
{
  biginteger mod = get_modulus(base, exp);
#ifdef DEBUG_bigmod
  if(mod.isNA() && !mpz_cmp_si(base.getValue().getValueTemp(), 1))
    Rprintf("bigmod pow(1, exp=%d)\n", mpz_get_si(exp.getValue().getValueTemp()));
  else if(mod.isNA() && !mpz_cmp_si(exp.getValue().getValueTemp(), 0))
    Rprintf("bigmod pow(base=%d, 0)\n", mpz_get_si(base.getValue().getValueTemp()));
#endif

  // if (base == 1  or  exp == 0)  return 1
  if(mod.isNA() &&
     ((!base.getValue().isNA() && !mpz_cmp_si(base.getValue().getValueTemp(), 1)) ||
      (! exp.getValue().isNA() && !mpz_cmp_si( exp.getValue().getValueTemp(), 0))))
    return DefaultBigMod(biginteger(1));
  if (base.getValue().isNA() || exp.getValue().isNA())
    return DefaultBigMod();
  int sgn_exp = mpz_sgn(exp.getValue().getValueTemp());
  bool neg_exp = (sgn_exp < 0); // b ^ -|e| =  1 / b^|e|
  mpz_t val; mpz_init(val); mpz_t_sentry val_s(val);
#ifdef DEBUG_bigmod
  Rprintf("bigmod pow(base=%3s, exp=%3s [mod=%3s]) ..\n",
          base.getValue().str(10).c_str(), exp.getValue().str(10).c_str(),
	  mod.str(10).c_str());
#endif
  if (mod.isNA()) { // <==> (both have no mod || both have mod. but differing)
    if(neg_exp) error(_("** internal error (negative powers for Z/nZ), please report!"));
    if (!mpz_fits_ulong_p(exp.getValue().getValueTemp()))
      error(_("exponent e too large for pow(z,e) = z^e"));// FIXME? return( "Inf" )
    // else :
    mpz_pow_ui(val, base.getValue().getValueTemp(),
 	       mpz_get_ui(exp.getValue().getValueTemp()));
  }
  else if( mpz_sgn(mod.getValueTemp()) != 0) { // check modulus non-zero
    if(neg_exp) { // negative exponent -- only ok if inverse exists
      if (mpz_invert(val, base.getValue().getValueTemp(), mod.getValueTemp()) == 0) {
	SEXP wOpt = Rf_GetOption1(Rf_install("gmp:warnNoInv"));
	if(wOpt != R_NilValue && Rf_asInteger(wOpt))
	  warning(_("pow(x, -|n|) returning NA as x has no inverse wrt modulus"));
	return(DefaultBigMod()); // return NA; was
      } // else: val = x^(-1) already: ==> result =  val ^ |exp| =  val ^ (-exp) :
      // nExp := - exp
      mpz_t nExp; mpz_init(nExp); mpz_neg(nExp, exp.getValue().getValueTemp());
      mpz_powm(val, val, nExp, mod.getValueTemp());
    } else { // non-negative exponent
      mpz_powm(val, base.getValue().getValueTemp(), exp.getValue().getValueTemp(), mod.getValueTemp());
    }
  }
  return DefaultBigMod(val, mod);
}

DefaultBigMod inv(const bigmod& x, const bigmod& m)
{
  if (x.getValue().isNA() || m.getValue().isNA())
    return DefaultBigMod();
  SEXP wOpt = Rf_GetOption1(Rf_install("gmp:warnNoInv"));
  bool warnI = (wOpt != R_NilValue && Rf_asInteger(wOpt));
  if (mpz_sgn(m.getValue().getValueTemp()) == 0) {
    if(warnI) warning(_("inv(0) returning NA"));
    return DefaultBigMod();
  }
  biginteger mod = get_modulus(x, m);
  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
  if (mpz_invert(val, x.getValue().getValueTemp(), m.getValue().getValueTemp()) == 0) {
    if(warnI) warning(_("inv(x,m) returning NA as x has no inverse modulo m"));
    return(DefaultBigMod()); // return NA; was
  }
  return DefaultBigMod(val, mod);
}

// R  as.bigz() :
DefaultBigMod set_modulus(const bigmod& x, const bigmod& m)
{
  if (!m.getValue().isNA() && mpz_sgn(m.getValue().getValueTemp()) == 0)
    error(_("modulus 0 is invalid"));
  //    if (!m.getValue().isNA() && mpz_cmp(x.getValue().getValueTemp(),m.getValue().getValueTemp())>=0) {
  if (!m.getValue().isNA() ) {
    DefaultBigMod t(x%m);
    return DefaultBigMod(t.getValue(), m.getValue());
  } else
    return DefaultBigMod(x.getValue(), m.getValue());
}

DefaultBigMod gcd(const bigmod& lhs, const bigmod& rhs)
{
  return create_bigmod(lhs, rhs, mpz_gcd);
}

DefaultBigMod lcm(const bigmod& lhs, const bigmod& rhs)
{
  return create_bigmod(lhs, rhs, mpz_lcm);
}


// return the modulus to use for the two bigmods.
// NA if incompatible.
biginteger get_modulus(const bigmod& b1, const bigmod& b2)
{
  if (b1.getModulus().isNA()) // NA: means "no modulus" <==> R's is.null(modulus(.))
    return b2.getModulus(); // if b2 is NA too, the return is correct: NA
  else if (b2.getModulus().isNA())
    return b1.getModulus();
  else if (mpz_cmp(b1.getModulus().getValueTemp(), b2.getModulus().getValueTemp())) {
    SEXP wOpt = Rf_GetOption1(Rf_install("gmp:warnModMismatch"));
    if(wOpt != R_NilValue && Rf_asInteger(wOpt))
      warning(_("modulus mismatch in bigz.* arithmetic"));
    return biginteger(); // i.e. NA
  } else // equal
    return b1.getModulus();
}


//    typedef void (*gmp_binary)(mpz_t, const mpz_t, const mpz_t);




// Create a bigmod from a binary combination of two other bigmods
DefaultBigMod create_bigmod(const bigmod& lhs, const bigmod& rhs, gmp_binary f,
		     bool zeroRhsAllowed) {
  if (lhs.getValue().isNA() || rhs.getValue().isNA())
    return DefaultBigMod();
  if (!zeroRhsAllowed && mpz_sgn(rhs.getValue().getValueTemp()) == 0) {
    warning(_("returning NA  for (modulus) 0 in RHS"));
    return DefaultBigMod();
  }
  biginteger mod = get_modulus(lhs, rhs);
  mpz_t val;
  mpz_init(val);
  mpz_t_sentry val_s(val);
  f(val, lhs.getValue().getValueTemp(), rhs.getValue().getValueTemp());
  //--- val := f(lhs, rhs)
#ifdef DEBUG_bigmod
  bool iNA = biginteger(val).isNA();
  char* buf;
  if(iNA)
    buf = NULL;
  else {
    buf = new char[mpz_sizeinbase(val, 10)+2];
    // possible minus sign, size of number + '\0'
    mpz_get_str(buf, 10, val);
  }
  Rprintf("create_bigmod(lhs=%3s, rhs=%3s [mod=%3s]) = %s%s",
	  lhs.getValue().str(10).c_str(),
	  rhs.getValue().str(10).c_str(),
	  mod.str(10).c_str(),
	  (iNA)? "NA" : buf,
	  (mod.isNA())? "\n" : " {before 'mod'}");
#endif
  if (!mod.isNA()) {
    mpz_mod(val, val, mod.getValueTemp());
#ifdef DEBUG_bigmod
    if(biginteger(val).isNA())
      Rprintf(" -> val = NA\n");
    else {
      mpz_get_str(buf, 10, val);
      Rprintf(" -> val = %s\n", buf);
    }
#endif
  }
  return DefaultBigMod(val, mod);
}
