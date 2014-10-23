#include <Rcpp.h>
#include "bbd.h"
using namespace Rcpp;

using RealType = double;

// [[Rcpp::export]]
std::vector< std::complex<double> > phi_Cpp (std::complex<double> s, int a0, int b0, std::vector<double> lambda2, std::vector<double> mu2, std::vector<double> x, std::vector<double> y, int A, int B) {  
  
  std::vector< std::complex<double> > phi((A+1-a0)*(B+1)*(B+1));
  const int dim = B+1+400;
  std::complex<double> fac,B1,B2,v;
  const std::complex<double> one(1.0,0.0), zero(0,0);
  std::vector<double> xvec(dim), prod_mu2((B+1)*(B+1)), prod_lambda2((B+1)*(B+1));
  std::vector< std::complex<double> > yvec(dim), lentz(B+1), Bk1dBk(B+1), BidBj((B+1)*(B+1));
        
      for (int a=0; a<=(A-a0); a++) {        
        for (int i=0; i<dim; i++) {
            xvec[i] = x[a+i*(A-a0+1)];
            yvec[i] = s+ y[a+i*(A-a0+1)];
        }
        
        lentz = lentz_Cpp(B,xvec,yvec);
        Bk1dBk = Bk1dBk_Cpp(B,xvec,yvec);
        BidBj = BidBj_Cpp(B,xvec,yvec,Bk1dBk);
        prod_mu2 = prod_vec_Cpp(a-a0+1,B,mu2);
        prod_lambda2 = prod_vec_Cpp(a-a0+1,B,lambda2);
        
        for (int i=0; i<=B; i++) 
        	for (int j=0; j<=B; j++) {
			      if (i<=j) {
				      if (i==j) {
					      fac = one;
				      } else {
					      fac = prod_mu2[(i+1)*(B+1)+j];
				      } 
				      B1 = BidBj[i*(B+1)+j];
				      B2 = one/Bk1dBk[j];
				      v = fac*B1/(B2+lentz[j]);
			      } else {
				      fac = prod_lambda2[j*(B+1)+i-1];
				      B1 = BidBj[j*(B+1)+i];
				      B2 = one/Bk1dBk[i];
				      v = fac*B1/(B2+lentz[i]);
			      }   
      phi[a+i*(A-a0+1)+j*(A-a0+1)*(B+1)] = v;
		}
  }
  return(phi);
}


