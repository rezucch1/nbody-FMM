/**
 * @file multipole_set3D.cpp
 * @brief Implementation of 3D MultipoleSet methods, 3D M2M translation, and 3D M2L conversion.
 */

#include "multipole_set3D.hpp"
#include "irregular.hpp"

constexpr std::complex<double> imaginary_power[] = {
    {1, 0},
    {0, 1},
    {-1, 0},
    {0, -1}
};

inline MultipoleSet<3>::MultipoleSet(unsigned int L)
: MultipoleSetI(L){
    elements.reserve(2*L*L + 4*L + 1);
}

std::complex<double> MultipoleSet<3>::operator()(unsigned l, int m) const{
    unsigned int idx = 2*l*l + 2*l + 2*m + 1;
    return std::complex<double>(elements[idx], elements[idx + 1]);
}

MultipoleSetI &MultipoleSet<3>::operator+=(const PowerSetI &z){ 
    elements[0] = z(0).real();
    unsigned int idx = 1;
    for (unsigned int n = 1; n <= L; ++n){
        for (int m = -n; m <= n; ++n){
            auto z_complex = z(n, m);
            elements[idx] = z_complex.real();
            elements[idx + 1] = z_complex.imag();
            idx += 2;
        }
    }
    return *(this);
}

std::unique_ptr<MultipoleSetI> MultipoleSet<3>::weigh_children_with_distance(const Tensor &d) const{

    auto M = new MultipoleSet<3>(L);
    M = 0;
    const auto regular = PowerSet<3>(L, d);

    M->elements[0] = elements[0];

    for (unsigned int j = 1; j <=L; ++j){
        for (int k = -j; k <= j; ++k){
            double A_row = 1;
            std::complex<double> m_complex = (*this)(0, 0) * A_row * regular(j, k);
            int m_row = k;
            for (unsigned int n = j-1; n >= 0; --n){
                if (m_row < 0){
                    ++m_row;
                    A_row *= std::sqrt(((2*n - j + k + 1)*(2*n - j + k + 2))/((2*j - 2*n)*(2*j - 2*n - 1)));
                } else if (m_row > 0){
                    --m_row;
                    A_row = std::sqrt(((2*n - j + k + 1)*(2*n - j + k + 2))/((2*j - 2*n)*(2*j - 2*n - 1)));
                } else { // m_row == 0
                    A_row = n + 1 / std::sqrt((j - k - n)*(j + k - n));
                }
                double A = A_row;
                m_complex += (*this)(j - n, k - m_row) * imaginary_power[(4 + (std::abs(k) - std::abs(m_row) - std::abs(k - m_row)) % 4) % 4] * A * regular(n, m_row);
                for (int m = m_row - 1; m >= -n && m >= k + n - j; --m){
                    A *= sqrt(((n + m + 1)*(j - k - n + m + 1)/(n - m)*(j + k - n - m)));
                    m_complex += (*this)(j - n, k - m) * imaginary_power[(4 + (std::abs(k) - std::abs(m) - std::abs(k - m)) % 4) % 4] * A * regular(n, m);
                }
                A = A_row;
                for (int m = m_row + 1; m <= n && m <= k + j - n; ++m ){
                    A *= sqrt(((n - m + 1)*(j + k - n - m + 1)/(n + m)*(j - k - n + m)));
                    m_complex += (*this)(j - n, k - m) * imaginary_power[(4 + (std::abs(k) - std::abs(m) - std::abs(k - m)) % 4) % 4] * A * regular(n, m);
                }
            }
            const unsigned int idx = 2*j*j + 2*j + 2*k + 1;
            M->elements[idx] = m_complex.real();
            M->elements[idx + 1] = m_complex.imag();
        }
    }
    return std::unique_ptr<MultipoleSetI>(M);
}

std::unique_ptr<LocalSetI> MultipoleSet<3U>::to_local(const Tensor &d) const{

    auto L = new LocalSet<3>(this->L);
    L = 0;
    const auto irregular = Irregular(2 * this->L, d);

    L->set_elements(0, 0, {elements[0], 0});

    for (unsigned int j = 1; j <= this->L; ++j){
        for (int k = -j; k <= j; ++k){
            double A_row = 1;
            int sign = 1;
            std::complex<double> l_complex = (*this)(0, 0) * A_row * irregular(j, k);
            for (unsigned int n = 1; n <= this->L; ++n){
                A_row = std::sqrt((n + j + k)*(n + j - k)) / n;
                sign = -sign;
                double A = A_row;
                l_complex += (double)sign * (*this)(n, 0) * A * irregular(j + n, k);
                for (int m = - 1; m >= -n; --m){
                    A *= sqrt(((n - m + j + k)*(n + m + 1)/(n + m + j - k + 1)*(n - m)));
                    l_complex += (double)sign * (*this)(n, m) * imaginary_power[(4 + (std::abs(k - m) - std::abs(k) - std::abs(m)) % 4) % 4] * A * irregular(j+n, k-m);
                }
                A = A_row;
                for (int m = 1; m <= n; ++m ){
                    A *= sqrt(((n + m + j - k)*(n - m + 1)/(n - m + j + k + 1)*(n + m)));
                    l_complex += (double)sign * (*this)(n, m) * imaginary_power[(4 + (std::abs(k - m) - std::abs(k) - std::abs(m)) % 4) % 4] * A * irregular(j+n, k-m);
                }
            }
            L->set_elements(j, k, l_complex);
        }
    }
    return std::unique_ptr<LocalSetI>(L);
}
