/**
 * @file multipole_set2D.cpp
 * @brief Implementation of 2D MultipoleSet M2M shifts and M2L conversion to local expansions.
 */

#include "multipole_set2D.hpp"

MultipoleSet<2>::MultipoleSet(unsigned int L)
: MultipoleSetI(L)
{
    elements.reserve(2*L + 1);
}

std::complex<double> MultipoleSet<2>::operator()(unsigned l, int m) const
{
    if (l == 0) return std::complex<double>(elements[0], 0);
    return std::complex<double>(elements[2*l-1], elements[2*l]);
}

std::unique_ptr<MultipoleSetI> MultipoleSet<2>::weigh_children_with_distance(const Tensor &d) const
{
    std::unique_ptr<MultipoleSet<2>> M = std::make_unique<MultipoleSet<2>>(L);
    *M = 0;
    M->elements[0] = MultipoleSetI::getElements(*this)[0];
    std::complex<double> d_k = 1;

    for(int k = 1; k <= L; k++){
        d_k *= std::complex{d[0], d[1]};
        unsigned int bin = 1;
        for( int l = k; l <= L; l++){
            std::complex mul_dk_cm = d_k * (*this)(l-k);
            M->elements[2*l-1] += bin * mul_dk_cm.real();
            M->elements[2*l] += bin * mul_dk_cm.imag();
            bin *= (l+1);
            bin /= (l+1-k);
        }
    }
    return M;
}

std::unique_ptr<LocalSetI> MultipoleSet<2>::to_local(const Tensor &d) const{
    
    std::unique_ptr<LocalSet<2>> local = std::make_unique<LocalSet<2>>(L);

    std::complex<double> z(d[0], d[1]);

    local->get_element(0) = (*this)(0)*std::log(z);

    std::complex<double> z_k = 1;

    for (unsigned int k = 1; k <= L; ++k){
        z_k *= z; 
        local->get_element(0) -= (*this)(k)/((double)k * z_k);
    }

    int sign = -1; 
    std::complex<double> z_n = 1;

    for (unsigned int n = 1; n <= L; ++n){
        sign = -sign;
        z_n *= z;
        local->get_element(n) = (double)sign * (*this)(0) / ((double) n * z_n);

        std::complex<double> z_kn = z_n;
        unsigned int bin_coef = 1;

        for (unsigned int k = 1; k <= L; ++k){
            z_kn *= z; 
            local->get_element(n) += (double)sign * bin_coef * (*this)(k)/((double)k * z_kn);
            bin_coef *= k + n;
            bin_coef /= k;
        }
    }

    return local;
}

template class MultipoleSet<2>;
