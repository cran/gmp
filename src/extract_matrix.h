/*! \file extract_matrix.h
 *  \brief functions to extract parts of matrix 
 *
 *  \version 1
 *
 *  \date Created: 25/06/11   
 *  \date Last modified: Time-stamp: <2006-06-17 23:08:24 antoine>
 *
 *  \author A. Lucas
 *
 * \note
 *  as usually, matrix x[i,j] (n x p) is represented by a vector
 *              x[i + j x n]  (i=0..n-1 ; j=0..p-1) 
 *
 *  \note Licence: GPL
 */


#ifndef EXTRACT_MATRIX_HEADER_GMP_R_
#define EXTRACT_MATRIX_HEADER_GMP_R_ 1



#include <R.h>
#include <Rdefines.h>

#include "bigvec_q.h"
#include "bigintegerR.h"


extern "C"
{

  /** \brief get subsets of a matrix */
  SEXP matrix_get_at_z(SEXP A,SEXP INDI, SEXP INDJ);

  /** \brief set subsets of a matrix */
  SEXP matrix_set_at_z(SEXP A,SEXP VAL ,SEXP INDI, SEXP INDJ);

  /** \brief get subsets of a matrix */
  SEXP matrix_get_at_q(SEXP A, SEXP INDI, SEXP INDJ);

  /** \brief set subsets of a matrix */
  SEXP matrix_set_at_q(SEXP A,SEXP VAL, SEXP INDI, SEXP INDJ);
 
  
}


namespace extract_gmp_R
{

  /** \brief Change R indices (in function like A[IND]) from
   * R (it can be logical/positive/negative) to a vector of size n
   * containing boolean: true: we update line/row, flase : we do not
   * change anything
   *
   */
  std::vector<bool> indice_set_at (unsigned int n , SEXP & IND);


  
  /**
   * \brief tranform a matrix from bigvec or bigvec_q format to 
   * vector<vector< big > > (a vector of lines)
   *
   * \note for bigvec: it does not take into account modulus. 
   *       
   */
  template< class T> std::vector<T> toVecVec(T& A)
  {
    std::vector<T> retour;
    unsigned int i;

    // case: vector
    if(A.nrow < 1)
      A.nrow = A.size();
    
    // check that size is a multiple of row
    if((A.size() / A.nrow) != static_cast<float>(A.size()) / static_cast<float>(A.nrow))
      Rf_error("malformed matrix");
    
    retour.resize(A.size() / A.nrow);
    for(i = 0; i < retour.size();  ++i)
      retour[i].value.resize(A.nrow);

    // go !
    for(i= 0 ; i < A.value.size(); ++i)
      // retour[col        ]  [row        ] 
      (retour[i / A.nrow ]).value[ i % A.nrow].setValue(A.value[i]);
	
    return(retour);

  }





  /**
   * \brief extract parts of a matrix
   *
   * \param A matrix (class bigvec or bigvec_q)
   * \param INDI,INDJ indices: either "LOGICAL" (true/false) or
   *  numbers: 
   *      - positives: we return row/col in INDI/INDJ
   *      - negatives: we retun all except row/col in INDI/INJ
   */
  template< class T> T get_at (T & A ,SEXP  & INDI, SEXP & INDJ)
  {
    // result = A[indii,indj]
    unsigned int oldnrow = A.nrow;

    std::vector<int> vi = bigintegerR::create_int(INDJ);

    std::vector<T> matricetmp = toVecVec(A);
    std::vector<T> matricetmp2;

    // only pointers
    std::vector<T> * matrice = &matricetmp;
    std::vector<T> * matricework = &matricetmp2; 

    T retour;

    unsigned int i,j,newnrow= 0 ;
    std::vector<int>::iterator it;

    // --------------------------
    // PART 1:  COLUMNS
    // --------------------------

    if(A.size()==0)
      return(retour);


    if(TYPEOF(INDJ) != NILSXP)	
      //LOCICAL
      if (TYPEOF(INDJ) == LGLSXP) 
	{
	    
	  // for all rows
	  unsigned int delta=0;
	  for (i = 0; i < (*matrice)[0].size(); ++i)
	    {
	      if (! vi[i+delta% vi.size()]) 
		{
		  matrice->erase(i+ matrice->begin());
		  --i; // indice in modified matrix
		  ++delta; // i+delta: old indice
		}
	    }
	  
	}
      else //INDJ: numbers
	{
	  vi.erase(remove(vi.begin(), vi.end(), 0),vi.end()) ; // remove all zeroes
	  if (vi.empty())
	    return retour;
	  
	  // case: a[-b]
	  // negatives...
	  if(vi[0] < 0)
	    {
	      // sort in reverse order
	      std::sort(vi.begin(),vi.end(),std::greater<int>() );
	      
	      // we should have vi like -1 -3 -7 -7 -12 ...

	      // remove duplicates
	      it = std::unique(vi.begin(),vi.end());
	      vi.erase(it,vi.end());

	      if ( vi.back() > 0)
		Rf_error("only 0's may mix with negative subscripts");


	      
	      it=vi.begin();
	      unsigned int delta=0;
	      // for all columns
	      for (j = 0; j < matrice->size(); ++j)
		{
		  if(it == vi.end())
		    break;

		  if (*it == - static_cast<int>(j+1+delta) ) 
		    {
		      matrice->erase(j+ matrice->begin());
		      ++it;
		      ++delta;
		      --j;
		    }
		  
		}
	      
	    }
	  else
	    // case: positive values: 1;5;7;7;9;10...
	    {

	      // note : we could have duplicate and indices are not
	      // sorted	
	      
	      
	      // allocate new matrix (all will be copied)
	      // number of columns
	      matricework->reserve(vi.size());

	      // for all [new] rows
	      for( it=vi.begin(); it != vi.end(); it++)
		{
		  if (*it  < 0)
		    Rf_error("only 0's may mix with negative subscripts");
		  if( static_cast<unsigned int>(*it-1) < matrice->size() )
		    {
		      //printf("on sort %s",(*matrice)[(*it)-1][0].str(10).c_str());
		      matricework->push_back( (*matrice)[*it-1] );
		    }
		      
		}
	      
	      // change addresses
	      matrice = &matricetmp2;
	      matricework = &matricetmp;
	      
	    }//end of case: int+positive values
	  
	}


    if(matrice->size()==0)
      return(retour);

    // --------------------------
    // PART 2:  ROWS
    // --------------------------
    vi = bigintegerR::create_int(INDI);
    if(TYPEOF(INDI) != NILSXP)	
      //LOCICAL
      if (TYPEOF(INDI) == LGLSXP) 
	{
	  // for all rows
	  unsigned int delta = 0;
	  for (i = 0; i < (*matrice)[0].size(); ++i)
	    {
	      if (! vi[(i+delta)% vi.size()]) 
		{
		  // for all columns
		  for (j = 0; j < matrice->size(); ++j)	     
		    (*matrice)[j].value.erase(i+(*matrice)[j].value.begin());
		      
		  //++newnrow;
		  --i; // i: new indice in modified matrix
		  ++delta; // i+delta = old indices
		}

	    }
	  
	}
      else //INDI: numbers
	{
	  vi.erase(remove(vi.begin(), vi.end(), 0),vi.end()) ; // remove all zeroes
	  if (vi.empty())
	    return retour;
	  
	  // case: a[-b]
	  // negatives...
	  if(vi[0] < 0)
	    {
	      std::sort(vi.begin(),vi.end(),std::greater<int>() );
	      // we should have vi like -1 -3 -7 -7 -12 ...

	      // remove duplicates
	      it = std::unique(vi.begin(),vi.end());
	      vi.erase(it,vi.end());
	    
	      if ( vi.back() > 0)
		Rf_error("only 0's may mix with negative subscripts");


	     
	      //newnrow = A.nrow;
	      it=vi.begin();
	      // for all rows
	      unsigned int delta = 0;
	      for (i = 0; i < (*matrice)[0].size(); ++i)
		{
		
		  if(it != vi.end() )
		    if (*it == - static_cast<int>(i+1+delta) ) 
		      {
			// for all columns
			for (j = 0; j < matrice->size(); ++j)
			  {
			    (*matrice)[j].value.erase(i+(*matrice)[j].value.begin());
			  }
			//--newnrow;
			--i; // i: new indice in modified matrix
			++delta; // i+delta = old indices
			++it;
		      }
		  
		}
	      
	    }
	  else
	    // case: positive values: 1;5;7;7;9;10...
	    {
	      // delete too large values
	      for(it = vi.begin(); it != vi.end(); ++it)
		if(*it > static_cast<int>((*matrice)[0].size()))
		  {
		    vi.erase(it);
		    --it;
		  }
		
	      // note : we could have duplicate and indices are not
	      // sorted	
	      
	      //newnrow = vi.size();
	      
	      // allocate new matrix (all will be copied)
	      // number of columns
	      matricework->resize( matrice->size());
	      // number of row
	      for (j = 0; j < matricework->size(); ++j)
		(*matricework)[j].resize( vi.size() );
	      

	      // for all [new] rows
	      for( i=0; i < vi.size(); ++i)
		{
		  if (vi[i]  < 0)
		    Rf_error("only 0's may mix with negative subscripts");
		  if( static_cast<unsigned int>(vi[i]-1) < (*matrice)[j].size() )
		    {		  
		      // for all columns
		      for (j = 0; j < matricework->size(); ++j)
			//newmat[i,j] = oldmat[ind[i],j]
			( (*matricework)[j]).value[i] = ((*matrice)[j]).value[vi[i]-1];
		    }
		  else
		    for (j = 0; j < matricework->size(); ++j)
		      ( (*matricework)[j]).value[i].setValue();
		}
	      
	      matrice = matricework; // change addresses
	    }//end of case: int+positive values
	  
	}


    // --------------------------
    // PART 3:  COPY RESULT
    // --------------------------

    newnrow =  (*matrice)[0].size();
    retour.resize(matrice->size() * newnrow);
    retour.nrow =  newnrow;
    for(i=0; i< newnrow ; ++i)
      for(j=0; j <  matrice->size() ; ++j)
	retour.value[i + j* newnrow ]  =  ((*matrice)[j]).value[i] ;
	
    if(oldnrow <= 0)
	retour.nrow = 0;

    return(retour);

  }  // end get_at


  /** \brief set a matrix: for R function src[idx,jdx] <- value
   *
   */
  template<class T> void set_at(T & src ,const T & value, SEXP & IDX, SEXP & JDX)
  {
    // case: vector
    if(src.nrow < 1)
      src.nrow = src.size();

    // check that size is a multiple of row
    if((src.size() / src.nrow) != static_cast<float>(src.size()) / static_cast<float>(src.nrow))
      Rf_error("malformed matrix");
    
    unsigned int ncol = src.size() / src.nrow; // number of col
    std::vector<bool> vidx =  indice_set_at ( src.nrow, IDX);
    std::vector<bool> vjdx =  indice_set_at ( ncol, JDX);
    
    unsigned int k=0;
    for(unsigned int j = 0 ; j < ncol; ++j)
      if(vjdx[j])
	{
	  for(unsigned int i = 0 ; i < src.nrow; ++i)
	    {
	      if(vidx[i])
		src.set(i + j * src.nrow, value[k % value.size()] );
	      ++k;
	    }
	}
    
    return;

  }//end set_at

}// end namespace






#endif
