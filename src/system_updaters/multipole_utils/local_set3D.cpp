/**
 * @file local_set3D.cpp
 * @brief Implementation of 3D LocalSet methods, 3D L2L translation, and 3D potential gradient evaluation.
 */

#include "local_set3D.hpp"
#include "power_set3D.hpp"
#include <cassert>

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
  assert(-(int)l <= m && m <= (int)l);
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
  assert(-(int)n <= m && m <= (int)n);
  unsigned int idx = 2*n*n + 2*n + 2*m;
  elements[idx] = l.real();
  elements[idx + 1] = l.imag();
}

std::unique_ptr<LocalSetI> LocalSet<3U>::distribute_parent_with_distance(const Tensor &d) const{
  auto L = std::make_unique<LocalSet<3>>(this->L);
  // *L = 0;
  const auto regular = PowerSet<3>(2 * this->L, d);

  for (unsigned int j = 0; j <= this->L; ++j){
    for (int k = -j; k <= (int)j; ++k){
      double A_row = 1;
      int sign = 1;
      std::complex<double> l_complex = (*this)(j, k) * regular(0, 0);
      for (unsigned int n = j + 1; n <= this->L; ++n){
        A_row *= std::sqrt((n - k)*(n + k)) / (n - j);
        sign = (n+j)%2 ? -1 : 1;
        double A = A_row;
        l_complex += (double)sign * (*this)(n, k) * A * regular(n - j, 0);
        for (int m = k - 1; m >= (int)(k - n + j); --m){
          A *= sqrt((double)((n - m)*(n + m -j - k + 1))/((n + m + 1)*(n - m - j + k)));
          l_complex += (double)sign * (*this)(n, m) * imaginary_power[(4 + (std::abs(m) - std::abs(m - k) - std::abs(k)) % 4) % 4] * A * regular(n - j, m - k);
        }
        A = A_row;
        for (int m = k + 1; m <= (int)(k + n - j); ++m ){
          A *= sqrt((double)((n + m)*(n - m - j + k + 1))/((n - m + 1)*(n + m - j - k)));
          l_complex += (double)sign * (*this)(n, m) * imaginary_power[(4 + (std::abs(m) - std::abs(m - k) - std::abs(k)) % 4) % 4] * A * regular(n - j, m - k);
        }
      }
      L->set_elements(j, k, l_complex);
      assert((*L)(j, k) == l_complex);
    }
  }
  return L;
}

Tensor LocalSet<3U>::get_gradient(const Tensor &d) const{
    Tensor gradient{0.0, 0.0, 0.0};

    PowerSet<3> regular(L, d);

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





/*Tensor LocalSet<3U>::get_gradient(const Tensor &d) const {

    auto evaluate_psi = [this](const Tensor &point) -> double {
        PowerSet<3> regular(this->L, point);
        std::complex<double> psi = 0.0;
        for (unsigned int n = 0; n <= this->L; ++n)
            for (int m = -static_cast<int>(n); m <= static_cast<int>(n); ++m)
                psi += (*this)(n, m) * regular(n, m);
        return psi.real();
    };

    const double eps = 1e-6 * std::max(1.0, d.norm());

    Tensor d_px = d; d_px[0] += eps;
    Tensor d_mx = d; d_mx[0] -= eps;
    Tensor d_py = d; d_py[1] += eps;
    Tensor d_my = d; d_my[1] -= eps;
    Tensor d_pz = d; d_pz[2] += eps;
    Tensor d_mz = d; d_mz[2] -= eps;

    Tensor gradient{0.0, 0.0, 0.0};
    gradient[0] = (evaluate_psi(d_px) - evaluate_psi(d_mx)) / (2*eps);
    gradient[1] = (evaluate_psi(d_py) - evaluate_psi(d_my)) / (2*eps);
    gradient[2] = (evaluate_psi(d_pz) - evaluate_psi(d_mz)) / (2*eps);

    return -gradient;
}

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

        return -gradient;
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
*/
