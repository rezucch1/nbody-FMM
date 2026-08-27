/**
 * @file power_set2D.cpp
 * @brief Implementation of 2D PowerSet complex power evaluation.
 */

#include "power_set2D.hpp"
#include "math_utils/math_utils.hpp"

PowerSet<2>::PowerSet(unsigned int L, const Tensor &r)
: PowerSetI(L)
{
    elements.reserve(2*L+1);
    elements.push_back(1.0);
    elements.push_back(r[0]);
    elements.push_back(r[1]);
    for(size_t i=3; i<2*L+1; i=i+2){
        std::complex<double> complex_multipole = std::complex(elements[i-2], elements[i-1]) * std::complex(r[0], r[1]);
        elements.push_back(complex_multipole.real());
        elements.push_back(complex_multipole.imag());
    } 
}

std::complex<double> PowerSet<2>::operator()(unsigned l, int m) const
{
    if (l == 0) return std::complex<double>(elements[0], 0);
    return std::complex<double>(elements[2*l-1], elements[2*l]);
}

template class PowerSet<2>;
