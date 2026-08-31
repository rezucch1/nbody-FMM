/**
 * @file irregular.cpp
 * @brief Implementation of Irregular solid harmonics calculation.
 */

#include "irregular.hpp"
#include <cmath>
#include <cassert>

Irregular::Irregular(unsigned int L, const Tensor &r)
: PowerSetI(L)
{
    
    elements.reserve((L+1)*(L+1));
    double r_squarenorm = r.squared_norm();

    elements.push_back(1 / std::sqrt(r_squarenorm));
    elements.push_back(r[2] * elements[0] / r_squarenorm); // z
    elements.push_back(-1/M_SQRT2 * r[0] * elements[0] / r_squarenorm);
    elements.push_back(-1/M_SQRT2 * r[1] * elements[0] / r_squarenorm);

    std::complex<double> el;

    for (int l=2; l<=L; l++){
        elements.push_back((r[2]*(2*l-1) * (*this)(l - 1, 0).real() - (l - 1) * (*this)(l - 2, 0).real()) / (l * r_squarenorm));
        for (int m=1; m<l-1; m++){
            el = ((2*l - 1) * r[2] * (*this)(l-1, m) - std::sqrt((l-1)*(l-1) - m*m) * (*this)(l-2, m)) / (std::sqrt(l*l - m*m) * r_squarenorm);
            elements.push_back(el.real());
            elements.push_back(el.imag());
        }

        // case m = l-1
        el = r[2] * std::sqrt(2*l - 1) * (*this)(l-1, l-1) / r_squarenorm;
        elements.push_back(el.real());
        elements.push_back(el.imag());

        // case m = l
        el = -sqrt((double)(2*l - 1)/(2*l)) * std::complex{r[0], r[1]} * (*this)(l-1, l-1) / r_squarenorm;
        elements.push_back(el.real());
        elements.push_back(el.imag());       
    }
}

std::complex<double> Irregular::operator()(unsigned l, int m) const
{
    assert(-(int)l <= m && m <= (int)l);
    if (m == 0)
        return std::complex<double>(elements[l*l], 0);
    else if (m > 0)
        return std::complex<double>(elements[l*l + 2*m - 1], elements[l*l + 2*m]);
    else {
        if (m % 2 == 0)
            return std::complex<double>(elements[l*l - 2*m - 1], -elements[l*l - 2*m]);
        else
            return std::complex<double>(-elements[l*l - 2*m - 1], elements[l*l - 2*m]);
    }
}
