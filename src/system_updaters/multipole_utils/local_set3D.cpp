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


/*
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
        (((*this)(j, -k) * std::conj(y_base)).real() + ((*this)(j, k) * y_base).real()) / denominator,
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
*/





Tensor LocalSet<3U>::get_gradient(const Tensor &d) const {
    Tensor gradient{0.0, 0.0, 0.0};

    const double x = d[0];
    const double y = d[1];
    const double z = d[2];

    const double r2 = x * x + y * y + z * z;
    const double r = std::sqrt(r2);

    // Caso 1: r = 0 (Origine) - Restituisce +grad(Phi)
    if (r < 1e-12) {
        const std::complex<double> L_1_1 = (*this)(1, 1);
        const double L_1_0 = (*this)(1, 0).real();

        // Segni del GRADIENTE POSITIVO (opposti alla forza)
        gradient[0] =  std::sqrt(2.0) * L_1_1.real(); 
        gradient[1] = -std::sqrt(2.0) * L_1_1.imag(); 
        gradient[2] =  L_1_0;                         

        return gradient;
    }

    const double rho2 = x * x + y * y;
    const double rho = std::sqrt(rho2);

    double dPhi_dr = 0.0;
    double dPhi_dtheta = 0.0;
    double dPhi_dphi = 0.0;

    PowerSet<3> regular(L, d); // R_n^m(d)

    for (unsigned int n = 0; n <= L; ++n) {
        for (int m = -static_cast<int>(n); m <= static_cast<int>(n); ++m) {
            
            const std::complex<double> L_nm = (*this)(n, m);
            const std::complex<double> R_nm = regular(n, m);

            // 1. dPhi/dr
            if (n > 0) {
                dPhi_dr += (L_nm * R_nm).real() * (static_cast<double>(n) / r);
            }

            // 2. dPhi/dphi
            if (m != 0) {
                const std::complex<double> val = L_nm * R_nm * std::complex<double>(0.0, static_cast<double>(m));
                dPhi_dphi += val.real();
            }

            // 3. dPhi/dtheta
            if (rho > 1e-12) {
                double term_theta = (z / r) * static_cast<double>(n) * (L_nm * R_nm).real();
                
                if (n > 0 && std::abs(m) <= static_cast<int>(n - 1)) {
                    const std::complex<double> R_n_1_m = regular(n - 1, m);
                    const double factor = std::sqrt(static_cast<double>(n * n - m * m));
                    term_theta -= (L_nm * R_n_1_m).real() * factor;
                }
                
                dPhi_dtheta += term_theta / rho;
            }
        }
    }

    // Proiezione in coordinate cartesiane (+grad(Phi))
    if (rho < 1e-12) {
        // Asse Z: rho = 0
        gradient[0] = 0.0;
        gradient[1] = 0.0;
        gradient[2] = dPhi_dr * (z >= 0.0 ? 1.0 : -1.0);
    } else {
        gradient[0] = (x / r) * dPhi_dr + (x * z / (r * rho)) * dPhi_dtheta - (y / (r*rho)) * dPhi_dphi;
        gradient[1] = (y / r) * dPhi_dr + (y * z / (r * rho)) * dPhi_dtheta + (x / (r*rho)) * dPhi_dphi;
        gradient[2] = (z / r) * dPhi_dr - (rho / r) * dPhi_dtheta;
    }

    return gradient;
}