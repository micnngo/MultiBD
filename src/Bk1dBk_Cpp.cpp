#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector< std::complex<double> > Bk1dBk_Cpp(int B, std::vector<double> xvec, std::vector< std::complex<double> > yvec) {
  
  std::vector< std::complex<double> > res(B+1);
  const std::complex<double> one(1,0), zero(0,0), tiny(1e-16,0);
  std::complex<double> aj, bj, Dj = zero, Dj1 = zero;
    for (int j=0; j<=B; j++) {
        aj = xvec[j];
        bj = yvec[j];
        Dj = bj + aj*Dj1;
        if (Dj == zero) Dj = tiny;
        Dj1 = one/Dj;
        res[j] = Dj1;
    }
    return(res);
}