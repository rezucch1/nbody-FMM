#include "local_set.hpp"
#include "local_set2D.hpp"


inline LocalSet<2>::LocalSet(unsigned int L) : LocalSetI(L){
  elements.reserve(L);
}

std::complex<double> LocalSet<2>::operator()(unsigned l, int m) const{
  return elements[l];
}

std::complex<double> &LocalSet<2>::get_element(unsigned int i)
{
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

LocalSet<2> *LocalSet<2>::distribute_parent_with_distance(const Tensor &d) const{
  std::complex<double> z(d[0], d[1]);
  auto child_local = new LocalSet<2>(L);

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

template class LocalSet<2>;