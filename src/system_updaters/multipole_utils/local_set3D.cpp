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
  unsigned int size = 2*L*L + 4*L + 2;
  elements.reserve(size);
  elements.assign(size, 0.0);
}

std::unique_ptr<LocalSetI> LocalSet<3U>::clone() const {
  return std::make_unique<LocalSet<3>>(*this);
}

std::complex<double> LocalSet<3U>::operator()(unsigned int l, int m) const{
  unsigned int idx = 2*l*l + 2*l + 2*m;
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
  unsigned int idx = 2*n*n + 2*n + 2*m;
  elements[idx] = l.real();
  elements[idx + 1] = l.imag();
}

std::unique_ptr<LocalSetI> LocalSet<3U>::distribute_parent_with_distance(const Tensor &d) const{
  auto L = std::make_unique<LocalSet<3>>(this->L);
  // *L = 0;
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
          A *= sqrt((double)((n - m)*(n + m -j - k + 1))/((n + m + 1)*(n - m - j + k)));
          l_complex += (double)sign * (*this)(n, m) * imaginary_power[(4 + (std::abs(m) - std::abs(m - k) - std::abs(k)) % 4) % 4] * A * regular(n - j, m - k);
        }
        A = A_row;
        for (int m = k + 1; m <= k + n - j; ++m ){
          A *= sqrt((double)((n + m)*(n - m - j + k + 1))/((n - m + 1)*(n + m - j - k)));
          l_complex += (double)sign * (*this)(n, m) * imaginary_power[(4 + (std::abs(m) - std::abs(m - k) - std::abs(k)) % 4) % 4] * A * regular(n - j, m - k);
        }
      }
      L->set_elements(j, k, l_complex);
    }
  }
  return L;
}

Tensor LocalSet<3U>::get_gradient(const Tensor &d) const
{
    Tensor gradient{0.0, 0.0, 0.0};

    PowerSet<3> regular(L, d);

for (unsigned int j = 0; j <= L; ++j) {
    for (int k = -static_cast<int>(j);
         k <= static_cast<int>(j);
         ++k) {

        auto coeff = (*this)(j, k);

        if (std::abs(coeff) > 1e-14) {
            std::cout << "  L(" << j << "," << k << ") = "
                      << coeff << std::endl;
        }
    }
}
    for (unsigned int j = 1; j <= L; ++j) {

        for (int k = -static_cast<int>(j);
             k <= static_cast<int>(j);
             ++k) {

            const std::complex<double> L_jk = (*this)(j, k);

            // dR_j^k / dx
            std::complex<double> dRdx{0.0, 0.0};
            if (k >= 0) {
              // R_{j-1}^{k-1}
              if (std::abs(k - 1) <= static_cast<int>(j - 1)) {
                dRdx -=
                0.5
                * std::sqrt(static_cast<double>((j + k) * (j + k - 1)))
                * regular(j - 1, k - 1);
              }
              // R_{j-1}^{k+1}
              if (std::abs(k + 1) <= static_cast<int>(j - 1)) {
                dRdx +=
                0.5
                * std::sqrt(static_cast<double>((j - k) * (j - k - 1)))
                * regular(j - 1, k + 1);
              }
            }
            else {
              // k < 0
              const int m = -k;
              std::complex<double> dRdx_pos{0.0, 0.0};
              if (std::abs(m - 1) <= static_cast<int>(j - 1)) {
                dRdx_pos -=
                0.5
                * std::sqrt(static_cast<double>((j + m) * (j + m - 1)))
                * regular(j - 1, m - 1);
              }
              if (std::abs(m + 1) <= static_cast<int>(j - 1)) {
                dRdx_pos +=
                0.5
                * std::sqrt(static_cast<double>((j - m) * (j - m - 1)))
                * regular(j - 1, m + 1);
              }
              dRdx = ((m % 2) ? -1.0 : 1.0)* std::conj(dRdx_pos);
}
// dR_j^k / dy
std::complex<double> dRdy{0.0, 0.0};

if (k >= 0) {

    if (std::abs(k - 1) <= static_cast<int>(j - 1)) {
        dRdy +=
            std::complex<double>(0.0, -0.5)
            * std::sqrt(
                static_cast<double>(
                    (j + k) * (j + k - 1)
                )
            )
            * regular(j - 1, k - 1);
    }

    if (std::abs(k + 1) <= static_cast<int>(j - 1)) {
        dRdy +=
            std::complex<double>(0.0, -0.5)
            * std::sqrt(
                static_cast<double>(
                    (j - k) * (j - k - 1)
                )
            )
            * regular(j - 1, k + 1);
    }
}
else {

    const int m = -k;

    std::complex<double> dRdy_pos{0.0, 0.0};

    if (std::abs(m - 1) <= static_cast<int>(j - 1)) {
        dRdy_pos +=
            std::complex<double>(0.0, -0.5)
            * std::sqrt(
                static_cast<double>(
                    (j + m) * (j + m - 1)
                )
            )
            * regular(j - 1, m - 1);
    }

    if (std::abs(m + 1) <= static_cast<int>(j - 1)) {
        dRdy_pos +=
            std::complex<double>(0.0, -0.5)
            * std::sqrt(
                static_cast<double>(
                    (j - m) * (j - m - 1)
                )
            )
            * regular(j - 1, m + 1);
    }

    dRdy =
        ((m % 2) ? -1.0 : 1.0)
        * std::conj(dRdy_pos);
}

std::complex<double> dRdz{0.0, 0.0};

if (std::abs(k) <= static_cast<int>(j - 1)) {
    dRdz =
        std::sqrt(
            static_cast<double>(
                (j - std::abs(k)) *
                (j + std::abs(k))
            )
        )
        * regular(j - 1, k);
}
            gradient[0] += (L_jk * dRdx).real();
            gradient[1] += (L_jk * dRdy).real();
            gradient[2] += (L_jk * dRdz).real();
        }
    }
    return gradient;
}

template class LocalSet<3>;