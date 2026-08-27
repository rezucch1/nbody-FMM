/**
 * @file local_set3D.cpp
 * @brief Implementation of 3D LocalSet methods, 3D L2L translation, and 3D potential gradient evaluation.
 */

#include "local_set3D.hpp"
#include "power_set3D.hpp"

constexpr std::complex<double> imaginary_power[] = {
  {1, 0},
  {0, 1},
  {-1, 0},
  {0, -1}
};

inline LocalSet<3>::LocalSet(unsigned int L) : LocalSetI(L){
  elements.reserve(2*L*L + 4*L + 1);
}

std::complex<double> LocalSet<3U>::operator()(unsigned int l, int m) const{
  unsigned int idx = 2*l*l + 2*l + 2*m + 2;
  return std::complex<double>(elements[idx], elements[idx + 1]);
}

LocalSet<3> &LocalSet<3>::operator+=(const LocalSetI *other){
  return (*this) += *(const LocalSet<3> *)other;
}

LocalSet<3> &LocalSet<3>::operator+=(const LocalSet<3> other){
  for (unsigned int i = 0; i < elements.size(); ++i)
    elements[i] += other.elements[i];

  return *this;
}

void LocalSet<3>::set_elements(unsigned int n, int m, std::complex<double> l){
  unsigned int idx = 2*n*n + 2*n + 2*m + 2;
  elements[idx] = l.real();
  elements[idx + 1] = l.imag();
}

LocalSet<3U> *LocalSet<3U>::distribute_parent_with_distance(const Tensor &d) const{
  auto L = new LocalSet<3>(this->L);
  L = 0;
  const auto regular = PowerSet<3>(2 * this->L, d);

  for (unsigned int j = 0; j <= this->L; ++j){
    for (int k = -j; k <= j; ++k){
      double A_row = 1;
      int sign = 1;
      std::complex<double> l_complex = (*this)(j, k) * regular(0, 0);
      for (unsigned int n = j + 1; n <= this->L; ++n){
        A_row *= std::sqrt((n - k)*(n + k)) / (n - j);
        sign = -sign;
        double A = A_row;
        l_complex += (double)sign * (*this)(n, k) * A * regular(n - j, 0);
        for (int m = k - 1; m >= k - n + j; --m){
          A *= sqrt(((n - m)*(n + m -j - k + 1)/(n + m + 1)*(n - m - j + k)));
          l_complex += (double)sign * (*this)(n, m) * imaginary_power[(4 + (std::abs(m) - std::abs(m - k) - std::abs(k)) % 4) % 4] * A * regular(n - j, k - m);
        }
        A = A_row;
        for (int m = k + 1; m <= k + n - j; ++m ){
          A *= sqrt(((n + m)*(n - m - j + k + 1)/(n - m + 1)*(n + m - j - k)));
          l_complex += (double)sign * (*this)(n, m) * imaginary_power[(4 + (std::abs(m) - std::abs(m - k) - std::abs(k)) % 4) % 4] * A * regular(n - j, k - m);
        }
      }
      L->set_elements(j, k, l_complex);
    }
  }
  return L;
}

Tensor LocalSet<3U>::get_gradient(const Tensor &d) const{
  Tensor gradient{0.0, 0.0, 0.0};
  PowerSet<3> regular(L, d);
  double denominator = d[0]*d[0] + d[1]*d[1];
  for (unsigned int j = 1; j <= L; ++j){
    gradient += {
      (*this)(j, 0).real() * regular(j, 0).real() * j * d[0] / denominator,
      (*this)(j, 0).real() * regular(j, 0).real() * j * d[1] / denominator,
      0.0
    };
    for (int k = 1; k <= j; ++k){
      const std::complex<double> x_base = regular(j, k) * std::complex<double>(d[0] * j, -d[1] * k);
      const std::complex<double> y_base = regular(j, k) * std::complex<double>(d[1] * j, d[0] * k);
      gradient += {
        (((*this)(j, -k) * std::conj(x_base)).real() + ((*this)(j, k) * x_base).real()) / denominator,
        (((*this)(j, -k) * std::conj(y_base)).real() + ((*this)(j, k) * x_base).real()) / denominator,
        0.0
      };
    }
    for (int k = -j + 1; k <= j - 1; ++k){
      const double z_gradient = ((*this)(j, k) * std::sqrt(j*j - k*k) * regular(j - 1, k)).real();
      gradient += {
        - z_gradient * d[0]*d[2] / denominator,
        - z_gradient * d[1]*d[2] / denominator,
        z_gradient
      };
    }
  }
  return gradient;
}

template class LocalSet<3>;