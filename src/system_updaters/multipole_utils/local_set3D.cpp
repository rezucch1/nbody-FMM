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
      const unsigned int idx = 2*j*j + 2*j + 2*k + 2;
      L->set_elements(j, k, l_complex);
    }
  }
  return L;
}

template class LocalSet<3>;