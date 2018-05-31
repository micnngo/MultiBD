#include "bbd.h"

template <class ParallelizationScheme>
std::vector<double> derivatives_lt_invert_Cpp_impl(double t,
    const std::vector<double>& lambda1, const std::vector<double>& lambda2,
    const double alpha, const double beta, const double powI_inf,
    const long int S0, const long int I0,
    const int Ap1, const int Bp1,
    const int ord, const int direction,
    const int nblocks, const double tol,
    ParallelizationScheme& scheme) {

//  auto start = std::chrono::steady_clock::now();

  const double double_PI  = 3.141592653589793238463, AA = 20.0;
  const int matsize = Ap1*Bp1;
  int kmax = nblocks;

  std::vector<mytype::ComplexVector> ig;
  std::vector<double> res(matsize), yvec(matsize);

  const size_t size = scheme.private_size();

  for (int i=0; i<matsize; ++i) {
    yvec[i] = lambda1[i] + lambda2[i];
  }

/////////////////////////////////////////////////////////////
// The following code computes the inverse Laplace transform
// Algorithm from Abate and Whitt using a Riemann sum
// Levin tranform is used to accelerate the convergence
/////////////////////////////////////////////////////////////

  ig.resize(kmax);

  scheme.for_each( boost::make_counting_iterator(0), boost::make_counting_iterator(kmax),
    [&](int w) {
      mytype::ComplexNumber s(AA/(2*t),double_PI*(w+1)/t);
      ig[w].resize(Ap1*Bp1);
      derivatives_lt_Cpp(s,lambda1,lambda2,alpha,beta,powI_inf,S0,I0,Ap1,Bp1,ord,direction,yvec,ig[w]);
    });

  mytype::ComplexNumber s(AA/(2*t),0.0);
  mytype::ComplexVector psum0(matsize);
  derivatives_lt_Cpp(s,lambda1,lambda2,alpha,beta,powI_inf,S0,I0,Ap1,Bp1,ord,direction,yvec,psum0);

  std::for_each(boost::make_counting_iterator(0), boost::make_counting_iterator(matsize),
    [&](int i) {
      Levin levin(tol); // A struct for Levin transform
      double term = 1e16, sdiff = 1e16;
      int k = 1;
      double psum = real(psum0[i])/(2*t);
      double sk,sk1;
      while ((std::abs(sdiff) > 1e-16)||(std::abs(term)>1e-3)) {
        double sgn = (k%2 == 0) ? 1.0 : -1.0;
        term = sgn*real(ig[k-1][i])/t;
        psum += term;
        double omega = k*term;
        sk = levin.next(psum,omega,1.0);
        if (k>1) sdiff = sk - sk1;
        k++;
        sk1 = sk;
        if (k > kmax) {
          ig.resize(kmax+nblocks);

          scheme.for_each( boost::make_counting_iterator(0), boost::make_counting_iterator(nblocks),
            [&](int w) {
              mytype::ComplexNumber s(AA/(2*t),double_PI*(w+kmax+1)/t);
              ig[w+kmax].resize(matsize);
              derivatives_lt_Cpp(s,lambda1,lambda2,alpha,beta,powI_inf,S0,I0,Ap1,Bp1,ord,direction,yvec,ig[w+kmax]);
            });

          kmax += nblocks;
        }
      }

      res[i] = sk1*exp(AA/2);
    });

//  auto end = std::chrono::steady_clock::now();
//
//  using TimingUnits = std::chrono::microseconds;
//  Rcpp::Rcout << "Time: " << std::chrono::duration_cast<TimingUnits>(end - start).count() << std::endl;

  return(std::move(res));
}

// [[Rcpp::export]]
std::vector<double> derivatives_lt_invert_Cpp(double t, const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const double alpha, const double beta, const double powI_inf,
    const long int S0, const long int I0,
    const int Ap1, const int Bp1, const int ord,
    const int direction, const int nblocks,
    const double tol, const int computeMode, const int nThreads) {

    switch(computeMode) {  // Run-time selection on compute_mode
      case 1: {
        loops::C11Threads loopC11Threads(nThreads, nblocks);
        return derivatives_lt_invert_Cpp_impl(t, lambda1, lambda2,
                    alpha, beta, powI_inf, S0, I0,
                    Ap1, Bp1, ord, direction,
                    nblocks, tol, loopC11Threads);
      }

      case 2: {
        loops::C11ThreadPool loopC11ThreadPool(nThreads, nblocks);
        return derivatives_lt_invert_Cpp_impl(t, lambda1, lambda2,
                    alpha, beta, powI_inf, S0, I0,
                    Ap1, Bp1, ord, direction,
                    nblocks, tol, loopC11ThreadPool);
      }

      case 3: {
        loops::C11Async loopC11Async(nThreads, nblocks);
        return derivatives_lt_invert_Cpp_impl(t, lambda1, lambda2,
                    alpha, beta, powI_inf, S0, I0,
                    Ap1, Bp1, ord, direction,
                    nblocks, tol, loopC11Async);
      }

      // case 4: {
      //   loops::RcppThreads loopRcppThreads(nThreads, nblocks);
      //   return derivatives_lt_invert_Cpp_impl(t, lambda1, lambda2,
      //               alpha, beta, S0, I0,
      //               Ap1, Bp1, ord, direction,
      //               nblocks, tol, loopRcppThreads);
      // }

      default: {
        loops::STL loopSTL;
        return derivatives_lt_invert_Cpp_impl(t, lambda1, lambda2,
                    alpha, beta, powI_inf, S0, I0,
                    Ap1, Bp1, ord, direction,
                    nblocks, tol, loopSTL);
      }
    }
}
