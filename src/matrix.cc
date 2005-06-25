/************************************************************/
/*! \file matrix.cc
 *  \brief C++ function to add matrix support 
 *
 *  \version 1
 *
 *  \date Created: 19/02/06   
 *  \date Last modified: Time-stamp: <2006-06-17 23:10:44 antoine>
 *
 *  \author A. Lucas
 *
 * \note
 *  as usually, matrix x[i,j] (n x p) is represented by a vector
 *              x[i + j x n]  (i=0..n-1 ; j=0..p-1) 
 *
 *  \note Licence: GPL
 */

#include <gmp.h>

#include <R.h>
#include <Rdefines.h>
//#include <list>

#include "biginteger.h"
//#include "biginteger_operator.h"
#include <vector>


#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("main", String)
#else
#define _(String) (String)
#endif


using namespace std;

#include "bigintegerR.h"
#include "matrix.h"


// C++ side of R function matrix()
SEXP as_matrixz (SEXP x, SEXP nrR, SEXP ncR, SEXP byrowR, SEXP mod)
{

  int nr, nc, byrow,i;
  
  //  int sizemod = matrixz_utils::length_modulus(x);


  nc=INTEGER(ncR)[0];
  nr=INTEGER(nrR)[0];
  byrow=INTEGER(byrowR)[0];

  // get "bigz" vector, this makes all conversion int to bigz etc...
  bigvec mat = bigintegerR::create_bignum(x);
  int lendat =mat.value.size();
  int sizemod = mat.modulus.size();
  // when modulus specified 
  bigvec modulus = bigintegerR::create_bignum(mod);
  if(modulus.value.size()>0) // should be allways the case
      if(!modulus.value[0].isNA())
	{
	  mat.modulus.resize(modulus.size());
	  for (i =0 ; i < (int)modulus.size(); i++)
	    mat.modulus[i] = modulus.value[i];
	  sizemod = modulus.size();
	}
  
  // A copy from R project to get correct dimension
  // all warnings...
  //
  if (nr == NA_INTEGER) // This is < 0 
    error(_("matrix: invalid 'nrow' value (too large or NA)"));
  if (nr < 0)
    error(_("matrix: invalid 'nrow' value (< 0)"));
  if (nc < 0)
    error(_("matrix: invalid 'ncol' value (< 0)"));
  if (nc == NA_INTEGER)
    error(_("matrix: invalid 'ncol' value (too large or NA)"));
  
  if(lendat > 0 ) {
    if (lendat > 1 && (nr * nc) % lendat != 0) {
      if (((lendat > nr) && (lendat / nr) * nr != lendat) ||
	  ((lendat < nr) && (nr / lendat) * lendat != nr))
	warning(_("data length [%d] is not a sub-multiple or multiple of the number of rows [%d] in matrix"), lendat, nr);
      else if (((lendat > nc) && (lendat / nc) * nc != lendat) ||
	       ((lendat < nc) && (nc / lendat) * lendat != nc))
	warning(_("data length [%d] is not a sub-multiple or multiple of the number of columns [%d] in matrix"), lendat, nc);
    }
    else if ((lendat > 1) && (nr * nc == 0)){
      warning(_("data length exceeds size of matrix"));
    }
  }

  // update dimension parameters 
  if(nr == 1)
    nr = (int)ceil(lendat / (double) nc);
  if(nc == 1)
    nc = (int)ceil(lendat / (double)nr);
    
  // when we extend "x"  
  if(nc*nr > lendat)
    {
      mat.value.resize(nr*nc);
      for(i = lendat; i< nr*nc; i++)
	mat.value[i] = mat.value[i % lendat];      
    }
  mat.nrow = nr;  
  if(byrow)
    {
      bigvec mat2 = matrixz::bigint_transpose (mat, nc,nr) ;
      mat2.nrow = nr;  
      return( bigintegerR::create_SEXP (mat2));
    }

  return( bigintegerR::create_SEXP (mat));
}


/*
 * Transposition
 */
SEXP bigint_transposeR(SEXP x)
{
  int nr=0, nc;
  SEXP dimAttr,dimName;
  //  int sizemod = matrixz_utils::length_modulus(x);

  PROTECT(dimAttr);
  PROTECT(dimName = Rf_allocVector(STRSXP,1));    
  SET_STRING_ELT(dimName, 0, Rf_mkChar("nrow"));
  dimAttr = Rf_getAttrib(x, dimName);
  UNPROTECT(2);
  if (TYPEOF(dimAttr) == INTSXP) {
    nr = INTEGER(dimAttr)[0];
  }
  else
    Rf_error("argument must be a matrix of class bigz.");

  bigvec mat = bigintegerR::create_bignum(x);
  nc = (int) mat.size() / nr;
  bigvec mat2 = matrixz::bigint_transpose(mat,nr,nc);
  //  mat2.nrow = nc;
  return( bigintegerR::create_SEXP( mat2));

}




/* \brief  matrix multiplication
 *  
 * returns matrix multiplication axb (a %*% b with R notation). 
 * \param a is of dimension nxp
 * \param b is of dimension pxm
 */
SEXP matrix_mul_z (SEXP a, SEXP b)
{
  int n, p, m;
  int i,j,k;
  int useMod=0;
  int dima[2], dimb[2];

  //  int sizemod_a = matrixz_utils::length_modulus(a);
  //int sizemod_b = matrixz_utils::length_modulus(b);

  bigvec mat_a = bigintegerR::create_bignum(a);
  bigvec mat_b = bigintegerR::create_bignum(b);

  int sizemod_a = mat_a.modulus.size();
  int sizemod_b = mat_b.modulus.size();

  dima[0] = mat_a.nrow;
  dimb[0] = mat_b.nrow;

  // cas: 2 vectors, then matrix multiplication returns scalar product 
  // first vector is allways considered as line-vector
  // second vector is allways considered as column-vector
  if(dima[0] < 1)
      dima[0] = 1;

  if(dimb[0] <= 1)
    dimb[0] = mat_b.size(); 

  dima[1] = mat_a.size() / dima[0];    
  dimb[1] = mat_b.size() / dimb[0];    

  /* dimension does not match */
  if(dima[1] != dimb[0])
    Rf_error("Size of matrix does not match");      
  
  n=dima[0];
  p=dima[1];
  m=dimb[1];
  

  //  cout << "mat "<< n<<"x"<<m<<"\n";
  bigvec res(n*m);
  res.nrow=n;

  mpz_t common_modulus ;
  mpz_t tmp_value ;

  mpz_init(tmp_value);
  mpz_init(common_modulus);


  /* modulus when modulus are "global" (i.e. of size 1) and 
   * either are the same, or only one of a or b is specified
   */
  if(! ((sizemod_a > 1) | (sizemod_b > 1)) )
    {
      if( (sizemod_a == 1 )&(sizemod_b == 1) )
	{
	  mpz_set(common_modulus, mat_a.modulus[0].getValueTemp());
	  if(mpz_cmp (common_modulus, mat_b.modulus[0].getValueTemp()) == 0)
	    useMod=1;      
	}
      else
	{
	  if(sizemod_a == 1 )
	    {
	      mpz_set(common_modulus,mat_a[0].modulus.getValueTemp());
	      useMod =1;
	    }
	  if( sizemod_b == 1)
	    {
	      mpz_set(common_modulus, mat_b[0].modulus.getValueTemp());
	      useMod =1;
	    }
	}
    }
  bigmod tmp;

  // here the computation:
  for(i=0; i<n; i++)
    for(j=0; j<m; j++)
      {
	res.value[i+ j*n].setValue(0);
	for(k=0; k<p; k++)
	  {
	    //	 tmp = create_bigmod(mat_a [i+ k*n_c  ],mat_b[k + j*p_c], mpz_mul);
	    mpz_mul(tmp_value, mat_a.value [i+ k*n  ].getValueTemp(),
		    mat_b.value[k + j*p].getValueTemp());
	    mpz_add(tmp_value,tmp_value,res.value[i+ j*n ].getValueTemp());
	    if(useMod)
	      mpz_mod(tmp_value,tmp_value,common_modulus);
	    res.value[i + j*n].setValue(tmp_value);	    
	  }
      }
  
  if(useMod)
    res.modulus.push_back(biginteger(common_modulus));

  mpz_clear(tmp_value);
  mpz_clear(common_modulus);

  return( bigintegerR::create_SEXP (res));  


}


SEXP biginteger_rbind(SEXP args) 
{
  int i=0,j=0; 
  bigvec result;
  bigvec v;

  result = bigintegerR::create_bignum(VECTOR_ELT(args,0));
  if(result.nrow ==0)
    result.nrow == result.size();

  result = matrixz::bigint_transpose(result,result.nrow,result.size() / result.nrow) ;

  for(i =1; i<LENGTH(args);i++)
    {
      v = bigintegerR::create_bignum(VECTOR_ELT(args,i));
      if(v.nrow == 0 )
	v.nrow = v.size();
      v = matrixz::bigint_transpose(v,v.nrow,v.size() / v.nrow) ;

      for(j=0; j< (int)v.size() ; j++)
	result.push_back(v[j]);
      v.clear();
    }

  result = matrixz::bigint_transpose(result,result.nrow,result.size() / result.nrow) ;

  return bigintegerR::create_SEXP(result);
}



namespace matrixz
{
  bigvec bigint_transpose ( bigvec & mat,int nr,int nc)
  {
    int i,j;

    /* cas: square matrix */
    bigvec matbis (nr * nc);
    matbis.nrow = nc;

    /* we compute transpose */
    for(i =0; i<nr; i++)
      for(j =0; j<nc; j++)
	matbis.set(j+i*nc,mat[i+j*nr]);
    
    return matbis;
  }


  /* return dimension in dima */
  int checkDims(int  dima,int  dimb)
  {
    if((dima>0) && ( dimb >0) )
      {
	if(  dimb != dima  )
	  Rf_error("Size of matrix does not match");      
	
      }
    else
      {
	/* either a or b is a matrix */
	if(dima == -1)
	  return(dimb);
      }
    return(dima);
  }
}
