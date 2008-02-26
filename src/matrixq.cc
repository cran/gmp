/************************************************************/
/*! \file matrixq.cc
 *  \brief C++ function to add matrix support 
 *
 *  \version 1
 *
 *  \date Created: 19/02/06   
 *  \date Last modified: Time-stamp: <2008-02-17 19:54:11 antoine>
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

#include "bigrational.h"
#include "bigrationalR.h"
//#include "biginteger_operator.h"
#include <vector>


#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("main", String)
#else
#define _(String) (String)
#endif


using namespace std;
//using namespace matrixz;

//#include "bigintegerR.h"
#include "matrixq.h"




// function called by matrix.bigz()
SEXP as_matrixq (SEXP x, SEXP nrR, SEXP ncR, SEXP byrowR, SEXP den)
{

  int nr, nc, byrow,i;
  
  //  int sizemod = matrixz_utils::length_modulus(x);


  nc=INTEGER(ncR)[0];
  nr=INTEGER(nrR)[0];
  byrow=INTEGER(byrowR)[0];

  /* get "bigz" vector, this makes all conversion int to bigz etc...*/
  bigvec_q mat = bigrationalR::create_bignum(x);
  int lendat =mat.size();

  bigvec_q denominator = bigrationalR::create_bignum(den);
  if(denominator.value.size()>0) // should be allways the case
      if(!denominator.value[0].isNA())
	{	  
	  for (i =0 ; i < (int)mat.size(); i++)
	    if(mat.value[i].isNA() && (denominator.value[i%denominator.size()].sgn() != 0) )
	      mat.value[i].setDenValue (denominator.value[i%denominator.size()].getValueTemp()) ; 

	}
  
  /* A copy from R project to get correct dimension
   * all warnings...
   */
  if (nr == NA_INTEGER) /* This is < 0 */
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

  /* update dimension parameters */
  if(nr == 1)
    nr = (int)ceil(lendat / (double) nc);
  if(nc == 1)
    nc = (int)ceil(lendat / (double)nr);
    
  /* when we extend "x"  */
  if(nc*nr > lendat)
    {
      mat.value.resize(nr*nc);
      for(i = lendat; i< nr*nc; i++)
	mat.value[i] = mat.value[i % lendat];      
    }
  mat.nrow = nr;  
  if(byrow)
   {
      bigvec_q mat2 = matrixq::bigq_transpose (mat, nc,nr) ;
      mat2.nrow = nr;  
      return( bigrationalR::create_SEXP (mat2));
   }

  return( bigrationalR::create_SEXP (mat));
}


// function called by t(m) when m is a bigrational
SEXP bigq_transposeR(SEXP x)
{
  int nr=0, nc;
  SEXP dimName;


  PROTECT(dimName = Rf_allocVector(STRSXP,1) );
  SET_STRING_ELT(dimName, 0, Rf_mkChar("nrow"));
  UNPROTECT(1);
  SEXP dimAttr = Rf_getAttrib(x, dimName);
 
  if (TYPEOF(dimAttr) == INTSXP) {
    nr = INTEGER(dimAttr)[0];
  }
  else
    Rf_error("argument must be a matrix of class bigq.");

  bigvec_q mat = bigrationalR::create_bignum(x);
  nc = (int) mat.size() / nr;
  bigvec_q mat_transp = matrixq::bigq_transpose(mat,nr,nc);
  mat_transp.nrow = nc;


  return( bigrationalR::create_SEXP( mat_transp));

}




/** matrix multiplication
 *  \brief returns matrix multiplication axb (a %*% b with R notation). 
 * a is of dimension nxp
 * b is of dimension pxm
 */
SEXP matrix_mul_q (SEXP a, SEXP b)
{
  int n, p, m;
  int i,j,k;
  int dima[2], dimb[2];

  //  int sizemod_a = matrixz_utils::length_modulus(a);
  //int sizemod_b = matrixz_utils::length_modulus(b);

  bigvec_q mat_a = bigrationalR::create_bignum(a);
  bigvec_q mat_b = bigrationalR::create_bignum(b);

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
  bigvec_q res(n*m);
  res.nrow=n;

  mpq_t tmp_value,product ;

  mpq_init(product);
  mpq_init(tmp_value);

  bigrational tmp;

  // here the computation:
  for(i=0; i<n; i++)
    for(j=0; j<m; j++)
      {
	mpq_set_ui(product,0,1);
	bool isna = false;
	for(k=0; k<p; k++)
	  {
	    // xij = sum aik bkj
	    if( !(mat_a[i+k*n].isNA() ||mat_b[ k + j * p].isNA()))
	      {
		mpq_mul(tmp_value, mat_a[i+k*n].getValueTemp(),mat_b[ k + j * p].getValueTemp());
		mpq_add(product,tmp_value,product);
	      }
	    else
	      {
		isna = true;
		break;
	      }
	  }
	if(isna)
	  {
	    res.value[i + j*n].setValue(0);
	    res.value[i + j*n].NA(true);
	  }
	else
	  res.value[i + j*n].setValue(product);
	    
      }
  
  

  mpq_clear(product);
  mpq_clear(tmp_value);
  return( bigrationalR::create_SEXP (res));  


}


SEXP bigrational_rbind(SEXP args) 
{
  int i=0,j=0; 
  bigvec_q result;
  bigvec_q v;

  result = bigrationalR::create_bignum(VECTOR_ELT(args,0));
  if(result.nrow ==0)
    result.nrow = result.size();

  result = matrixq::bigq_transpose(result,result.nrow,result.size() / result.nrow) ;

  for(i =1; i<LENGTH(args);i++)
    {
      v = bigrationalR::create_bignum(VECTOR_ELT(args,i));
      if(v.nrow == 0 )
	v.nrow = v.size();
      v = matrixq::bigq_transpose(v,v.nrow,v.size() / v.nrow) ;

      for(j=0; j< (int)v.size() ; j++)
	result.push_back(v[j]);
      v.clear();
    }

  result = matrixq::bigq_transpose(result,result.nrow,result.size() / result.nrow) ;

  return bigrationalR::create_SEXP(result);
}


bigvec_q matrixq::bigq_transpose (const  bigvec_q & mat,int nr,int nc)
{
  int i,j;
  bigvec_q matbis (nr * nc);
  matbis.nrow = nc;

  /* we compute transpose */
  for(i =0; i<nr; i++)
    for(j =0; j<nc; j++)
	matbis.value[j+i*nc].setValue(mat.value[i+j*nr]);


  return(matbis);
}
 

