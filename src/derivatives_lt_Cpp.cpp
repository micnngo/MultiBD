#include "bbd.h"

void derivatives_lt_Cpp(const mytype::ComplexNumber s, const std::vector<double>& lambda1,
    const std::vector<double>& lambda2, const double alpha, const double beta,
    const long int S0, const long int I0, const int Ap1, const int Bp1,
    const int ord, const int direction,
    const std::vector<double>& yvec, mytype::ComplexVector& ff) {

/////////////////////////////////////////////////////////////////////////////////
// The following code computes Laplace transform of the transition probabilities
/////////////////////////////////////////////////////////////////////////////////

  mytype::ComplexVector f(Ap1*Bp1);

  if (direction == 0) { // Forward direction
    f[0] = 1/(s + yvec[0]);
    switch(ord) {
      case 1: {
        ff[0] = - I0*f[0]/(s + yvec[0]);
        break;
      }

      case 2: {
        ff[0] = - S0*I0*f[0]/(s + yvec[0]);
        break;
      }
    }

    for (int i=0; i<(Ap1-1); ++i) {
      f[(i+1)*Bp1] = lambda1[i]*f[i*Bp1]/(s + yvec[i+1]);
      switch(ord) {
        case 1: {
          ff[(i+1)*Bp1] = (lambda1[i]*ff[i*Bp1] - lambda2[i+1]*f[(i+1)*Bp1]/alpha)/(s + yvec[i+1]);
          break;
        }

        case 2: {
          ff[(i+1)*Bp1] = (lambda1[i]*ff[i*Bp1] + lambda1[i]*f[i*Bp1]/beta - lambda1[i+1]*f[(i+1)*Bp1]/beta)/(s + yvec[i+1]);
          break;
        }
      }
    }

    for (int j=0; j<(Bp1-1); ++j) {
      f[j+1] = lambda2[j*Ap1]*f[j]/(s + yvec[(j+1)*Ap1]);
      switch(ord) {
        case 1: {
          ff[j+1] = (lambda2[j*Ap1]*ff[j] + lambda2[j*Ap1]*f[j]/alpha - lambda2[(j+1)*Ap1]*f[j+1]/alpha)/(s + yvec[(j+1)*Ap1]);
          break;
        }

        case 2: {
          ff[j+1] = (lambda2[j*Ap1]*ff[j] - lambda1[(j+1)*Ap1]*f[j+1]/beta)/(s + yvec[(j+1)*Ap1]);
          break;
        }
      }
    }

    for (int i=0; i<(Ap1-1); ++i) {
      for (int j=0; j<(Bp1-1); ++j) {
        f[(i+1)*Bp1 + j+1] = (lambda1[i + (j+1)*Ap1]*f[i*Bp1 + j+1] + lambda2[i+1 + j*Ap1]*f[(i+1)*Bp1 + j])/(s + yvec[i+1 + (j+1)*Ap1]);
        switch(ord) {
          case 1: {
            ff[(i+1)*Bp1 + j+1] = (lambda1[i + (j+1)*Ap1]*ff[i*Bp1 + j+1] + lambda2[i+1 + j*Ap1]*ff[(i+1)*Bp1 + j] + lambda2[i+1 + j*Ap1]*f[(i+1)*Bp1 + j]/alpha - lambda2[i+1 + (j+1)*Ap1]*f[(i+1)*Bp1 + j+1]/alpha)/(s + yvec[i+1 + (j+1)*Ap1]);
            break;
          }

          case 2: {
            ff[(i+1)*Bp1 + j+1] = (lambda1[i + (j+1)*Ap1]*ff[i*Bp1 + j+1] + lambda2[i+1 + j*Ap1]*ff[(i+1)*Bp1 + j] + lambda1[i + (j+1)*Ap1]*f[i*Bp1 + j+1]/beta - lambda1[i+1 + (j+1)*Ap1]*f[(i+1)*Bp1 + j+1]/beta)/(s + yvec[i+1 + (j+1)*Ap1]);
            break;
          }
        }
      }
    }


  } else { // Backward direction
    f[Ap1*Bp1-1] = 1/(s + yvec[Ap1*Bp1-1]);

    for (int i=(Ap1-2); i>=0; --i) {
      f[i*Bp1 + Bp1-1] = lambda1[i + (Bp1-1)*Ap1]*f[(i+1)*Bp1 + Bp1-1]/(s + yvec[i+ (Bp1-1)*Ap1]);
    }

    for (int j=(Bp1-2); j>=0; --j) {
      f[(Ap1-1)*Bp1 + j] = lambda2[Ap1-1 + j*Ap1]*f[(Ap1-1)*Bp1 + j+1]/(s + yvec[Ap1-1 + j*Ap1]);
    }

    for (int i=(Ap1-2); i>=0; --i) {
      for (int j=(Bp1-2); j>=0; --j) {
        f[i*Bp1 + j] = (lambda1[i + j*Ap1]*f[(i+1)*Bp1 + j] + lambda2[i + j*Ap1]*f[i*Bp1 + j+1])/(s + yvec[i + j*Ap1]);
      }
    }
  }

}
