/************************************************************/
/*! \file gmpR.c
 *  \brief C function to interface R and libgmp
 *
 *  \version 1
 *
 *  \date 26/09/04   
 *  \date Last modified: Time-stamp: <2004-09-27 16:21:42 lucas>
 *
 *  \author A. Lucas
 *
 *  \note Licence: CeCIll
 */


#include <stdlib.h>
#include <math.h>
#include <gmp.h>

#include <R.h>
#include <Rdefines.h>


/* Globals variables */

gmp_randstate_t seed_state;
int seed_init=0;

/* Functions that return a boolean, integer or double*/ 


/** @brief Check if nbrein is prime, Use mpz_probab_prime_p from libgmp
           return 2 if prime, 0 if not, 1 if don't know
    @param nbrein Integer to check (as a character)
    @param val      Return value (0,1,2)
    @param reps     Number of repeat
*/
isprime (char** nbrein,int * val,int* reps)
{
  mpz_t    a;
  /*mpz_init (a);*/
  /* On copie le nombre dans le format de la lib gmp */
  mpz_init_set_str(a,*nbrein,10) ;

 
  *val = mpz_probab_prime_p(a,*reps);
  mpz_clear (a);

}

/** @brief division as a double. res = (doulbe) a / b
    @param a     Integer
    @param b     Integer
    @param res   Result: a/b
*/
int gmp_div (char ** a, char ** b, double * res)
{

  mpq_t    az;
  mpq_t    bz;
  mpq_init(az);
  mpq_init(bz);

  mpq_set_str(az,*a,10) ;
  mpq_set_str(bz,*b,10) ;

  /* Here: division */
  mpq_div(az,az,bz);
  
  /* Return double */
  *res = mpq_get_d(az);

  /* free memory */
  mpq_clear (az);
  mpq_clear (bz);


}


/* Functions of NxN -> N  */ 


/** @brief Give next prime number
    @param nbrein   Number input
    @param nbreout  Prime number (out)
*/
SEXP nextprime(SEXP nbrein,SEXP nbreout)
{
  mpz_t    a,b;
  int l;
  char * mot;
  /*  mpz_init (a);   */
  mpz_init (b);

  PROTECT (nbrein = AS_CHARACTER(nbrein));
  /* STRING_ELT(nbrein,0) est l'element 0 du "vecteur nbrein
     CHAR(..) est un char *
     printf("%s\n",CHAR(STRING_ELT(nbrein,0)));*/
  mpz_init_set_str(a,CHAR(STRING_ELT(nbrein,0)),10) ;
  UNPROTECT(1);

  mpz_nextprime(b,a);
  l = mpz_sizeinbase(b,10)+2;
  
  mot = (char *) malloc(l * sizeof(char));
  if(mot == NULL)
    {
      printf("GMP: Memory allocation error at nextprime\n");
      return (0);
    }
  mpz_get_str(mot,10,b);

  /* Ecriture dans le vecteur de R 
  PROTECT (nbreout = AS_CHARACTER(nbreout));
  SET_STRING_ELT (nbreout,0,mkChar(mot)); 
  UNPROTECT(1);*/

  PROTECT (nbreout = allocVector(STRSXP,1));
  SET_STRING_ELT (nbreout,0,mkChar(mot));  
  UNPROTECT(1);

  free(mot);
  
  mpz_clear (a);
  mpz_clear (b);

  return(nbreout);
}



/** @brief Addition of a & b
    @param a Integer
    @param b Integer
*/
SEXP gmp_add (SEXP a, SEXP b)
{
  SEXP ans;
  mpz_t    az;
  mpz_t    bz;
  char * mot;
  int l;

  /* store input data into appropriate mode */
  /*mpz_init (az);
    mpz_init (bz);*/
  PROTECT (a = AS_CHARACTER(a));
  PROTECT (b = AS_CHARACTER(b));
  /* STRING_ELT(a,0) est l'element 0 du "vecteur nbrein
     CHAR(..) est un char *
     printf("%s\n",CHAR(STRING_ELT(a,0)));*/
  mpz_init_set_str(az,CHAR(STRING_ELT(a,0)),10) ;
  mpz_init_set_str(bz,CHAR(STRING_ELT(b,0)),10) ;
  UNPROTECT(2);

  /* Here: addition */
  mpz_add(az,az,bz);
  
  l = mpz_sizeinbase(az,10)+2;  
  mot = (char *) malloc(l * sizeof(char));
  if(mot == NULL)
    {
      printf("GMP: Memory allocation error at gmp_add\n");
      return (0);
    }
  mpz_get_str(mot,10,az);
  
  PROTECT (ans = allocVector(STRSXP,1));
  SET_STRING_ELT (ans,0,mkChar(mot));  
  UNPROTECT(1);

  mpz_clear (az);
  mpz_clear (bz);
  free(mot);
  return(ans);

}


/** @brief Substraction of a & b
    @param a Integer
    @param b Integer
*/
SEXP gmp_sub (SEXP a, SEXP b)
{
  SEXP ans;
  mpz_t    az;
  mpz_t    bz;
  char * mot;
  int l;

  /* store input data into appropriate mode */
  /*mpz_init (az);
    mpz_init (bz);*/
  PROTECT (a = AS_CHARACTER(a));
  PROTECT (b = AS_CHARACTER(b));
  /* STRING_ELT(a,0) est l'element 0 du "vecteur nbrein
     CHAR(..) est un char *
     printf("%s\n",CHAR(STRING_ELT(a,0)));*/
  mpz_init_set_str(az,CHAR(STRING_ELT(a,0)),10) ;
  mpz_init_set_str(bz,CHAR(STRING_ELT(b,0)),10) ;
  UNPROTECT(2);

  /* Here: substraction */
  mpz_sub(az,az,bz);
  
  l = mpz_sizeinbase(az,10)+2;  
  mot = (char *) malloc(l * sizeof(char));
  if(mot == NULL)
    {
      printf("GMP: Memory allocation error at gmp_sub\n");
      return (0);
    }
  mpz_get_str(mot,10,az);
  
  PROTECT (ans = allocVector(STRSXP,1));
  SET_STRING_ELT (ans,0,mkChar(mot));  
  UNPROTECT(1);
  mpz_clear (az);
  mpz_clear (bz);
  free(mot);
  return(ans);

}


/** @brief Multiplication of a & b
    @param a Integer
    @param b Integer
*/
SEXP gmp_mul (SEXP a, SEXP b)
{
  SEXP ans;
  mpz_t    az;
  mpz_t    bz;
  char * mot;
  int l;

  /* store input data into appropriate mode */
  /*mpz_init (az);
    mpz_init (bz);*/
  PROTECT (a = AS_CHARACTER(a));
  PROTECT (b = AS_CHARACTER(b));
  /* STRING_ELT(a,0) est l'element 0 du "vecteur nbrein
     CHAR(..) est un char *
     printf("%s\n",CHAR(STRING_ELT(a,0)));*/
  mpz_init_set_str(az,CHAR(STRING_ELT(a,0)),10) ;
  mpz_init_set_str(bz,CHAR(STRING_ELT(b,0)),10) ;
  UNPROTECT(2);

  /* Here: multiplication */
  mpz_mul(az,az,bz);
  
  l = mpz_sizeinbase(az,10)+2;  
  mot = (char *) malloc(l * sizeof(char));
  if(mot == NULL)
    {
      printf("GMP: Memory allocation error at gmp_mul\n");
      return (0);
    }
  mpz_get_str(mot,10,az);
  
  PROTECT (ans = allocVector(STRSXP,1));
  SET_STRING_ELT (ans,0,mkChar(mot));  
  UNPROTECT(1);
  mpz_clear (az);
  mpz_clear (bz);
  free(mot);
  return(ans);

}

/** @brief Quotient of a / b
    @param a Integer
    @param b Integer
*/
SEXP gmp_fdivq (SEXP a, SEXP b)
{
  SEXP ans;
  mpz_t    az;
  mpz_t    bz;
  char * mot;
  int l;

  /* store input data into appropriate mode */
  /*  mpz_init (az);
      mpz_init (bz);*/
  PROTECT (a = AS_CHARACTER(a));
  PROTECT (b = AS_CHARACTER(b));
  /* STRING_ELT(a,0) est l'element 0 du "vecteur nbrein
     CHAR(..) est un char *
     printf("%s\n",CHAR(STRING_ELT(a,0)));*/
  mpz_init_set_str(az,CHAR(STRING_ELT(a,0)),10) ;
  mpz_init_set_str(bz,CHAR(STRING_ELT(b,0)),10) ;
  UNPROTECT(2);

  /* Here: Division q */
  mpz_fdiv_q(az,az,bz);
  
  l = mpz_sizeinbase(az,10)+2;  
  mot = (char *) malloc(l * sizeof(char));
  if(mot == NULL)
    {
      printf("GMP: Memory allocation error at gmp_fdivq\n");
      return (0);
    }
  mpz_get_str(mot,10,az);
  
  PROTECT (ans = allocVector(STRSXP,1));
  SET_STRING_ELT (ans,0,mkChar(mot));  
  UNPROTECT(1);
  mpz_clear (az);
  mpz_clear (bz);
  free(mot);
  return(ans);

}

/** @brief Reminder of a / b (a = bq +r)
    @param a Integer
    @param b Integer
*/
SEXP gmp_fdivr (SEXP a, SEXP b)
{
  SEXP ans;
  mpz_t    az;
  mpz_t    bz;
  char * mot;
  int l;

  /* store input data into appropriate mode */
  /*  mpz_init (az);
      mpz_init (bz);*/
  PROTECT (a = AS_CHARACTER(a));
  PROTECT (b = AS_CHARACTER(b));
  /* STRING_ELT(a,0) est l'element 0 du "vecteur nbrein
     CHAR(..) est un char *
     printf("%s\n",CHAR(STRING_ELT(a,0)));*/
  mpz_init_set_str(az,CHAR(STRING_ELT(a,0)),10) ;
  mpz_init_set_str(bz,CHAR(STRING_ELT(b,0)),10) ;
  UNPROTECT(2);

  /* Here: substraction */
  mpz_fdiv_r(az,az,bz);
  
  l = mpz_sizeinbase(az,10)+2;  
  mot = (char *) malloc(l * sizeof(char));
  if(mot == NULL)
    {
      printf("GMP: Memory allocation error at gmp_fdivr\n");
      return (0);
    }
  mpz_get_str(mot,10,az);
  
  PROTECT (ans = allocVector(STRSXP,1));
  SET_STRING_ELT (ans,0,mkChar(mot));  
  UNPROTECT(1);
  mpz_clear (az);
  mpz_clear (bz);
  free(mot);
  return(ans);

}


/** @brief Great common divisor of a & b
    @param a Integer
    @param b Integer
*/
SEXP gmp_gcd (SEXP a, SEXP b)
{
  SEXP ans;
  mpz_t    az;
  mpz_t    bz;
  char * mot;
  int l;

  /* store input data into appropriate mode */
  PROTECT (a = AS_CHARACTER(a));
  PROTECT (b = AS_CHARACTER(b));
  /* STRING_ELT(a,0) est l'element 0 du "vecteur nbrein
     CHAR(..) est un char *
     printf("%s\n",CHAR(STRING_ELT(a,0)));*/
  mpz_init_set_str(az,CHAR(STRING_ELT(a,0)),10) ;
  mpz_init_set_str(bz,CHAR(STRING_ELT(b,0)),10) ;
  UNPROTECT(2);

  /* Here: substraction */
  mpz_gcd(az,az,bz);
  
  l = mpz_sizeinbase(az,10)+2;  
  mot = (char *) malloc(l * sizeof(char));
  if(mot == NULL)
    {
      printf("GMP: Memory allocation error at gmp_gcd\n");
      return (0);
    }
  mpz_get_str(mot,10,az);
  
  PROTECT (ans = allocVector(STRSXP,1));
  SET_STRING_ELT (ans,0,mkChar(mot));  
  UNPROTECT(1);
  mpz_clear (az);
  mpz_clear (bz);
  free(mot);
  return(ans);

}



/** @brief Lower common multiple of a & b
    @param a Integer
    @param b Integer
*/
SEXP gmp_lcm (SEXP a, SEXP b)
{
  SEXP ans;
  mpz_t    az;
  mpz_t    bz;
  char * mot;
  int l;

  /* store input data into appropriate mode */
  PROTECT (a = AS_CHARACTER(a));
  PROTECT (b = AS_CHARACTER(b));
  /* STRING_ELT(a,0) est l'element 0 du "vecteur nbrein
     CHAR(..) est un char *
     printf("%s\n",CHAR(STRING_ELT(a,0)));*/
  mpz_init_set_str(az,CHAR(STRING_ELT(a,0)),10) ;
  mpz_init_set_str(bz,CHAR(STRING_ELT(b,0)),10) ;
  UNPROTECT(2);

  /* Here: substraction */
  mpz_lcm(az,az,bz);
  
  l = mpz_sizeinbase(az,10)+2;  
  mot = (char *) malloc(l * sizeof(char));
  if(mot == NULL)
    {
      printf("GMP: Memory allocation error at gmp_lcm\n");
      return (0);
    }
  mpz_get_str(mot,10,az);
  
  PROTECT (ans = allocVector(STRSXP,1));
  SET_STRING_ELT (ans,0,mkChar(mot));  
  UNPROTECT(1);
  mpz_clear (az);
  mpz_clear (bz);
  free(mot);
  return(ans);

}


/** @brief Invert of a in [b/bZ]
    @param a Integer
    @param b Integer
*/
SEXP gmp_invert (SEXP a, SEXP b)
{
  SEXP ans;
  mpz_t    az;
  mpz_t    bz;
  char * mot;
  int l,res;

  /* store input data into appropriate mode */
  PROTECT (a = AS_CHARACTER(a));
  PROTECT (b = AS_CHARACTER(b));
  /* STRING_ELT(a,0) est l'element 0 du "vecteur nbrein
     CHAR(..) est un char *
     printf("%s\n",CHAR(STRING_ELT(a,0)));*/
  mpz_init_set_str(az,CHAR(STRING_ELT(a,0)),10) ;
  mpz_init_set_str(bz,CHAR(STRING_ELT(b,0)),10) ;
  UNPROTECT(2);

  /* Here: Calculus */
  res = mpz_invert(az,az,bz);
  
  if(res==0)
    {
      mot = (char *) malloc(2 * sizeof(char));
      mot[0]='0';
      mot[1]='\0';
    }
  else
    {
      l = mpz_sizeinbase(az,10)+2;  
      mot = (char *) malloc(l * sizeof(char));
      if(mot == NULL)
	{
	  printf("GMP: Memory allocation error at gmp_lcm\n");
	  return (0);
	}
      mpz_get_str(mot,10,az);
    }

  PROTECT (ans = allocVector(STRSXP,1));
  SET_STRING_ELT (ans,0,mkChar(mot));  
  UNPROTECT(1);
  mpz_clear (az);
  mpz_clear (bz);
  free(mot);
  return(ans);

}


/** @brief Power of a ^ exp [m]
    @param a Integer
    @param exp Integer
    @param m Integer
*/
SEXP gmp_powm (SEXP a, SEXP exp, SEXP m)
{
  /* ans = a ^ exp [ m ] */
  SEXP ans;
  mpz_t    az;
  mpz_t    ez;
  mpz_t    mz;
  char * mot;
  int l;

  /* store input data into appropriate mode */
  PROTECT (a = AS_CHARACTER(a));
  PROTECT (exp = AS_CHARACTER(exp));
  PROTECT (m = AS_CHARACTER(m));
  /* STRING_ELT(a,0) est l'element 0 du "vecteur nbrein
     CHAR(..) est un char *
     printf("%s\n",CHAR(STRING_ELT(a,0)));*/
  mpz_init_set_str(az,CHAR(STRING_ELT(a,0)),10) ;
  mpz_init_set_str(ez,CHAR(STRING_ELT(exp,0)),10) ;
  mpz_init_set_str(mz,CHAR(STRING_ELT(m,0)),10) ;
  UNPROTECT(3);

  /* Here: calculus */
  mpz_powm(az,az,ez,mz);
  
  l = mpz_sizeinbase(az,10)+2;  
  mot = (char *) malloc(l * sizeof(char));
  if(mot == NULL)
    {
      printf("GMP: Memory allocation error at gmp_lcm\n");
      return (0);
    }
  mpz_get_str(mot,10,az);
  
  PROTECT (ans = allocVector(STRSXP,1));
  SET_STRING_ELT (ans,0,mkChar(mot));  
  UNPROTECT(1);
  mpz_clear (az);
  mpz_clear (ez);
  mpz_clear (mz);
  free(mot);
  return(ans);

}


/** @brief Power of a ^ exp
    @param a   Integer
    @param exp Integer
*/
SEXP gmp_pow (SEXP a, SEXP exp)
{
  /* ans = a ^ exp  */
  SEXP ans;
  mpz_t    az;
  mpz_t    ez;
  char * mot;
  int l;
  long int e ;

  /* store input data into appropriate mode */
  PROTECT (a = AS_CHARACTER(a));
  PROTECT (exp = AS_INTEGER(exp));
  /* STRING_ELT(a,0) est l'element 0 du "vecteur nbrein
     CHAR(..) est un char *
     printf("%s\n",CHAR(STRING_ELT(a,0)));*/
  mpz_init_set_str(az,CHAR(STRING_ELT(a,0)),10) ;
  e = INTEGER(exp)[0];

  UNPROTECT(2);

  /* Here: calculus */
  mpz_pow_ui(az,az,e);
  
  l = mpz_sizeinbase(az,10)+2;  
  mot = (char *) malloc(l * sizeof(char));
  if(mot == NULL)
    {
      printf("GMP: Memory allocation error at gmp_lcm\n");
      return (0);
    }
  mpz_get_str(mot,10,az);
  
  PROTECT (ans = allocVector(STRSXP,1));
  SET_STRING_ELT (ans,0,mkChar(mot));  
  UNPROTECT(1);
  mpz_clear (az);

  free(mot);
  return(ans);

}

/** @brief Bezoult coefficients: compute g,s and t as as + bt = g 
    @param a Integer
    @param b Integer
*/
SEXP gmp_gcdex (SEXP a, SEXP b)
{
  SEXP ans;
  mpz_t    az;
  mpz_t    bz;
  mpz_t    gz;
  mpz_t    sz;
  mpz_t    tz;
  char * mot1;
  char * mot2;
  char * mot3;
  int l1,l2,l3;


  /* store input data into appropriate mode */
  PROTECT (a = AS_CHARACTER(a));
  PROTECT (b = AS_CHARACTER(b));
  /* STRING_ELT(a,0) est l'element 0 du "vecteur nbrein
     CHAR(..) est un char *
     printf("%s\n",CHAR(STRING_ELT(a,0)));*/
  mpz_init_set_str(az,CHAR(STRING_ELT(a,0)),10) ;
  mpz_init_set_str(bz,CHAR(STRING_ELT(b,0)),10) ;
  UNPROTECT(2);

  mpz_init(gz);
  mpz_init(sz);
  mpz_init(tz);
  /* Here: calculus */
  mpz_gcdext(gz,sz,tz,az,bz);
  mpz_clear(az);
  mpz_clear(bz);


  l1 = mpz_sizeinbase(gz,10)+2;  
  mot1 = (char *) malloc(l1 * sizeof(char));
  if(mot1 == NULL)
    {
      printf("GMP: Memory allocation error at gmp_lcm\n");
      return (0);
    }
  mpz_get_str(mot1,10,gz);
  mpz_clear(gz);

  l2 = mpz_sizeinbase(sz,10)+2;  
  mot2 = (char *) malloc(l2 * sizeof(char));
  if(mot2 == NULL)
    {
      printf("GMP: Memory allocation error at gmp_lcm\n");
      return (0);
    }
  mpz_get_str(mot2,10,sz);
  mpz_clear(sz);

  l3 = mpz_sizeinbase(tz,10)+2;  
  mot3 = (char *) malloc(l3 * sizeof(char));
  if(mot3 == NULL)
    {
      printf("GMP: Memory allocation error at gmp_lcm\n");
      return (0);
    }
  mpz_get_str(mot3,10,tz);
  mpz_clear(tz);
  

  PROTECT (ans = allocVector(STRSXP,3));
  SET_STRING_ELT (ans,0,mkChar(mot1));  
  SET_STRING_ELT (ans,1,mkChar(mot2));  
  SET_STRING_ELT (ans,2,mkChar(mot3));  
  UNPROTECT(1);
  

  free(mot1);
  free(mot2);
  free(mot3);
  return(ans);

}

/** @brief Random number generation
    \note If seed is not initialised: generation of a new seed
    @param length Integer number will be of length 2^length
    @param newseed Integer, seed initialisation (if exists)
    @param ok Integer 1: seed generation 0 not
*/
SEXP gmp_rand_u (SEXP length,SEXP newseed, SEXP ok)
{
  SEXP ans;
  mpz_t    az;
  mpz_t    bz;
  char * mot;
  int l,flag,len;
  
  extern int seed_init;
  extern gmp_randstate_t seed_state;


  /* store input data into appropriate mode */
  PROTECT (newseed = AS_CHARACTER(newseed));
  PROTECT (ok = AS_INTEGER(ok));
  PROTECT (length = AS_INTEGER(length));
  flag = INTEGER(ok)[0];
  len = INTEGER(length)[0];
  if(flag == 1)
    mpz_init_set_str(az,CHAR(STRING_ELT(newseed,0)),10) ;

  UNPROTECT(3);

  /*  printf("seed, %d\n",seed_init);*/
  /* Random seed initialisation */

  if(seed_init==0)
    {
      gmp_randinit_default(seed_state);
      printf("Seed default initialisation\n");
    }
  if(flag == 1)
    {
      gmp_randseed(seed_state,az);
      printf("Seed initialisation\n");
    }

  seed_init = 1;


  /*  Random number generation  */
  mpz_init (bz);
  mpz_urandomb(bz,seed_state,len);
  
  l = mpz_sizeinbase(bz,10)+2;  
  mot = (char *) malloc(l * sizeof(char));
  if(mot == NULL)
    {
      printf("GMP: Memory allocation error at gmp_lcm\n");
      return (0);
    }
  mpz_get_str(mot,10,bz);
  
  PROTECT (ans = allocVector(STRSXP,1));
  SET_STRING_ELT (ans,0,mkChar(mot));  
  UNPROTECT(1);

  mpz_clear (bz);
  if (flag == 1)
    mpz_clear (az);

  free(mot);
  return(ans);

}

