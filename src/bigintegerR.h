

namespace bigintegerR{

  vector<bigmod> create_vector(SEXP param);

  vector<bigmod> create_bignum(SEXP param);

  SEXP create_SEXP(const vector<bigmod>& v);


}
