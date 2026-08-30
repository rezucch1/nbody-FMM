/**
 * @file local_set2D.cpp
 * @brief Implementation of 2D LocalSet methods, L2L translation, and 2D potential gradient evaluation.
 */

#include "local_set.hpp"
#include "local_set2D.hpp"

inline LocalSet<2>::LocalSet(unsigned int L) : LocalSetI(L){
  elements.reserve(L + 1);
  elements.assign(L + 1, {0.0, 0.0});
}

std::unique_ptr<LocalSetI> LocalSet<2>::clone() const {
  return std::make_unique<LocalSet<2>>(*this);
}

std::complex<double> LocalSet<2>::operator()(unsigned l, int m) const{
  return elements[l];
}

std::complex<double> &LocalSet<2>::get_element(unsigned int i)
{
  return elements[i];
}

const std::complex<double> &LocalSet<2>::get_element(unsigned int i) const{
  return elements[i];
}

LocalSet<2> &LocalSet<2>::operator+=(const LocalSetI *other){
  return (*this) += *(const LocalSet<2> *)other;
}

LocalSet<2> &LocalSet<2>::operator+=(const LocalSet<2> other){
  for (unsigned int l = 0; l <= L; ++l)
    elements[l] += other.elements[l];

  return *this;
}

std::unique_ptr<LocalSetI> LocalSet<2>::distribute_parent_with_distance(const Tensor &d) const{
  std::complex<double> z(d[0], d[1]);
  auto child_local = std::make_unique<LocalSet<2>>(L);

  for (unsigned int k = 0; k <= L; ++k){
    unsigned int bin_coef = 1;
    std::complex<double> z_nk = 1;
    for (unsigned int n = k; n <= L; ++n){
      child_local->elements[k] += (double)bin_coef * elements[n] * z_nk;
      z_nk *= z;
      bin_coef *= n + 1;
      bin_coef /= n + 1 - k;
    }
  }

  return child_local;
}

Tensor LocalSet<2U>::get_gradient(const Tensor &d) const{
  Tensor grad{0.0, 0.0};
  auto pow = std::complex(1.0, 0.0); 
  for (unsigned int n = 1; n <= L; ++n){
    std::complex<double> complex_grad_n = (double)n * get_element(n) * pow;
    grad += {-complex_grad_n.real(), +complex_grad_n.imag()};
    pow *= std::complex<double>{d[0], d[1]};
  }
  return grad;
}

template class LocalSet<2>;