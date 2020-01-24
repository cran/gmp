
#include "bigvec.h"

/** \brief constructor
 *
 */
bigvec::bigvec(unsigned int i) :
  math::Matrix<bigmod>(),
  value(i),
  modulus(0),
  valuesMod(),
  nrow(-1)
{
}


bigvec::bigvec(const bigvec & vecteur) :
  math::Matrix<bigmod>(),
  value(vecteur.value.size()),
  modulus(vecteur.modulus.size()),
  valuesMod(),
  nrow(vecteur.nrow)
{
  //  *this = vecteur;
  value.resize(vecteur.value.size());
  modulus.resize(vecteur.modulus.size());
  std::vector<biginteger>::const_iterator jt=vecteur.modulus.begin();
  std::vector<biginteger>::iterator it = modulus.begin();
  while(it != modulus.end())
    {
      *it = *jt;
      ++it;
      ++jt;
    }
  jt = vecteur.value.begin();
  for(it=value.begin(); it != value.end(); ++it)
    {
      *it= *jt;
      ++jt;
    }
}

bigvec::~bigvec(){
  clearValuesMod();
}

//
std::string bigvec::str(int i,int b) const
{
    if (value[i].isNA())
      return ("NA");

    std::string s; // sstream seems to collide with libgmp :-(
    bool onemodulus = modulus.size()>0;
    if(onemodulus)
      onemodulus = !modulus[i%modulus.size()].isNA() ;
    if (onemodulus)
	s = "(";
    s += value[i].str(b);
    if (onemodulus) {
	s += " %% ";
	s += modulus[i%modulus.size()].str(b);
	s += ")";
    }

    return s;
}

bigmod & bigvec::get(unsigned int row, unsigned int col) {
  return (*this)[row + col*nrow];
}


 bigmod & bigvec::operator[] (unsigned int i) 
{
  checkValuesMod();
  return *valuesMod[i];
}


const bigmod & bigvec::operator[] (unsigned int i) const
{
  bigvec * nonConst = const_cast<bigvec*>(this);
  nonConst->checkValuesMod();
  return *valuesMod[i];
}

void bigvec::set(unsigned int row, unsigned int col, const  bigmod & val) {
  set( row + col*nrow,val);
}

void bigvec::checkValuesMod() {
  if (valuesMod.size() != value.size()){
    // reconstruct bigmod that are references to values and modulus:
    clearValuesMod();
    if(modulus.size()>0) {
      for (int i = 0 ; i < value.size(); i++)
	valuesMod.push_back(new bigmod(value[i], modulus[i%modulus.size()]));
    } else {
      for (int i= 0 ; i < value.size(); i++)
	valuesMod.push_back(new BigModInt(value[i]));
    }
    
  }
}

void bigvec::clearValuesMod(){
  for (int i = 0 ; i < valuesMod.size(); i++){
    delete valuesMod[i];
  }
  valuesMod.clear();
}

void bigvec::set(unsigned int i,const bigmod & val)
{
  value[i] = val.getValue();

  if(!val.getModulus().isNA())
    {
      if(modulus.size()> i)
	{
	  modulus[i] = val.getModulus() ;
	  return;
	}

      int nrow_mod = nrow;
      if(nrow<1)
	nrow_mod = 1;
      if( (modulus.size() ==  (unsigned int)nrow_mod ) || (modulus.size() == 1) )
	{
	  // check "row-mod" or "global" modulus
	  if(!(val.getModulus() != modulus[i % modulus.size()])) {
	    return;
	  }
	}

      // we add "previous"
      nrow_mod = modulus.size();
      for(unsigned int j = nrow_mod; j< i;++j)
	modulus.push_back(modulus[j % nrow_mod]);
      modulus.push_back(val.getModulus());
      clearValuesMod();
    }
}

void bigvec::push_back(const bigmod & number)
{
  int nrow_mod = (nrow < 0) ? 1 : nrow;
  clearValuesMod();

  value.push_back(number.getValue());

  if((!number.getModulus().isNA()) || (modulus.size()>0) )
    {
      // pathologic case: value.size >0 ; modulus.size =0
      // we assume previous mod where NA
      if((modulus.size() ==0) && (value.size()>0))
	{
	  modulus.resize(value.size()-1);
	  modulus.push_back( number.getModulus());
	  return;
	}

      // standard cas
      if((modulus.size() != 1 ) && (static_cast<int>(modulus.size()) != nrow_mod) )
	{
	  modulus.push_back(number.getModulus());
	  return;
	}

      // pathologic case:
      //  value modulus [nrow=4]
      //  2     2  push_back: ok
      //  2     2  push_back: nothing
      //  2     1  push_back: shoud add previous modulus then 1
      // note nrow_mod is either 1 ither nrow (when nrow >1)
      nrow_mod = modulus.size();
      if(modulus[(value.size() -1)% nrow_mod ] != number.getModulus())
	{
	  // we add "previous"
	  for(unsigned int i = nrow_mod; i< value.size()-1;i++)
	    modulus.push_back(modulus[i % nrow_mod]);
	  modulus.push_back(number.getModulus());
	  }
    }
}

/** 
 * insert int value
 */
void bigvec::push_back(int value_p)
{
  clearValuesMod();
  value.push_back(biginteger(value_p));
}

/** 
 * insert int value
 */
void bigvec::push_back(biginteger & value_p)
{
  clearValuesMod();
  value.push_back(value_p);
}

/**
 * Insert Big Integer value
 */
void bigvec::push_back(const __mpz_struct * value_p)
{
  clearValuesMod();
  value.push_back(biginteger(value_p));
}


// return size of value
unsigned int bigvec::size() const
{
  return(value.size());
}


unsigned int  bigvec::nRows() const {
   return abs(nrow);
}


// hummm. to avoid !
void bigvec::resize(unsigned int i)
{
  clearValuesMod();


  value.resize(i);
  if(i < modulus.size())
      modulus.resize(i);
}

// clear all
void bigvec::clear()
{
  clearValuesMod();
  value.clear();
  modulus.clear();
  nrow = -1;
}


// assignment operator
bigvec & bigvec::operator= (const bigvec & rhs)
{
  if(this != &rhs)
    {
      value.resize(rhs.value.size());
      modulus.resize(rhs.modulus.size());
      std::vector<biginteger>::const_iterator jt=rhs.modulus.begin();
      std::vector<biginteger>::iterator it = modulus.begin();
      while(it != modulus.end())
	{
	  *it = *jt;
	  ++it;
	  ++jt;
	}
      jt = rhs.value.begin();
      for(it=value.begin(); it != value.end(); ++it)
	{
	  *it= *jt;
	  ++jt;
	}
      nrow = rhs.nrow;
    }
  return(*this);
}



// Comparison operator
bool operator!=(const bigvec & rhs, const bigvec& lhs)
{
  std::vector<biginteger>::const_iterator it = rhs.value.begin();
  std::vector<biginteger>::const_iterator jt = lhs.value.begin();

  if( (rhs.value.size() != lhs.value.size()) || \
      (rhs.nrow != lhs.nrow )  )
    return(false);

  // check value
  while(it != rhs.value.end())
    {
      if(*it != *jt)
	return(false);
      it++; jt++;
    }
  for(unsigned int i=0; i < std::max( rhs.modulus.size() ,lhs.modulus.size() ); ++i)
    if(rhs.modulus[i % rhs.modulus.size()] != \
       lhs.modulus[i % lhs.modulus.size()] )
      return(false);
  return(true);
}



// never used
void bigvec::print()
{
  if(nrow > 0) {
    for(int i=0; i < nrow; ++i)
      {
	for(unsigned int j=0; j < (value.size() / nrow); ++j)
	  Rprintf("%s\t", value[i+j* nrow].str(10).c_str() );
	Rprintf("\n");
      }
  }
  else {
    for(unsigned int i=0; i < value.size(); ++i)
      Rprintf("%s\t", value[i].str(10).c_str() );
    Rprintf("\n");
  }
}
