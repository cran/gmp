#ifndef TEMPLATE_MATRIX_GMP
#define TEMPLATE_MATRIX_GMP 1



namespace matrix{

  template< class T>
    class Matrix{
    /*

    class Row{
    private :
      Matrix & source;
      int start;
      Row(Matrix & source_p, int rowIndex)
	: source(source_p),
	  start(rowIndex*source.nCols()){};
      ~Row(){};
    public:
      T  operator[] (unsigned int i) const {
	return source.get(start + i);
      };


    };
    */

  private:
    Matrix * transposate;

  public: 
    Matrix();

    virtual ~Matrix() {
      if (transposate) delete transposate;
    };
    /*
    Row  get (unsigned int i) const{
      Row row (*this, i);
      return row;
      }*/

    /** i, j -> index =  i + j*rows */
    //  virtual T  operator[](unsigned int index) const=0;

    /** return numRow*runCols */
    virtual unsigned int size() const = 0;

    virtual unsigned int nRows() const = 0;
    
    virtual unsigned int nCols() const{
      return size() / nRows();
    };


    virtual  T & get(unsigned int row, unsigned int col) = 0;

    virtual  void set(unsigned int row, unsigned int col, const  T& val) =0;

    /**
     * \brief assign a value at indice i
     */
    //  virtual void set(unsigned int i,const T & val) = 0;
    
    /**
     * \brief substract lambda * line j to line i
     */
    void subLine(unsigned int i,unsigned int j,const T & lambda){
      for(int k=0; k < nCols() ; ++k)
	set(i , k ,  get(i, k) - ( get(j , k) * lambda ) );
    }

    /**
     * \brief multiply line i by lambda
     */
    void mulLine(unsigned int i,const T & lambda){
      for(int k=0; k < nCols(); ++k)
	set(i , k,  get(i,k)  * lambda  );
    }
    
    Matrix & transpose();
  };
  
  template< class T>
    class Transpose: public Matrix<T> {
    private :
      Matrix<T> & source;
      
     
    public:
      Transpose(Matrix<T> & source_p)
	: source(source_p) {
      };

      unsigned int size(){
	return source.size();
      }

      unsigned int nRows(){
	return source.nCols();
      }

      T & get(unsigned int i, unsigned int j) const{
	return source.get(j,i);
      }

      void set(unsigned int i,unsigned int j, const T & val){
	source.set(j,i,val);
      }

    };

  // template code.
  template<class T>
    Matrix<T>::Matrix() : transposate(NULL)
  {
  }

  template<class T>
    Matrix<T> & Matrix<T>::transpose(){
    if (transposate == 0){
      transposate =new Transpose<T>(*this);
    }
    return *transposate;
  }
  
}

#endif
