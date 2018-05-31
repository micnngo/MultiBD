#include "bbd.h"

// [[Rcpp::export]]
std::vector<double> SIR_derivatives_Cpp(const double t, //const long int N, //don't needed?
    const double alpha, const double beta, const double powI_inf, const double powI_rem,
    const long int S0, const long int I0, const int Ap1, const int Bp1,
    const int ord, const int direction,
    const int nblocks, const double tol, const int computeMode, const int nThreads) {

      const int matsize = Ap1*Bp1;
      std::vector<double> lambda1(matsize), lambda2(matsize);

      for (int a=0; a<Ap1; ++a) {
        for (int b=0; b<Bp1; ++b) {
          double Spop = std::max(0.0, S0-a+0.0);
          double Ipop = std::max(0.0, a+I0-b+0.0);
          lambda1[a + b*Ap1] = beta*Spop*pow(Ipop, powI_inf); // Infection rate is beta*S*I
          lambda2[a + b*Ap1] = alpha*pow(Ipop, powI_rem);    //Recovery is gamma*I
        }
      }

      return(derivatives_lt_invert_Cpp(t, lambda1, lambda2,
            alpha, beta, powI_inf, powI_rem, S0, I0,
            Ap1, Bp1, ord, direction,
            nblocks, tol, computeMode, nThreads));
}
