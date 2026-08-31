/**
 * @file power_set3D.cpp
 * @brief Implementation of 3D PowerSet regular solid harmonics evaluation.
 */

#include "power_set3D.hpp"
#include <cmath>
#include <cassert>

PowerSet<3>::PowerSet(unsigned int L, const Tensor &r)
: PowerSetI(L)
{
    elements.reserve((L+1)*(L+1));
    elements.push_back(1);
    elements.push_back(r[2]); // z
    elements.push_back(-1/M_SQRT2 * r[0]);
    elements.push_back(-1/M_SQRT2 * r[1]);

    double r_squarenorm = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    std::complex<double> el;

    for (int l=2; l<=L; l++){
        elements.push_back((r[2]*(2*l-1) * (*this)(l - 1, 0).real() - (l - 1) * r_squarenorm * (*this)(l - 2, 0).real()) / (double)l);
        assert((*this)(l, 0).imag() == 0);
        for (int m=1; m<l-1; m++){
            el = ((2*l - 1) * r[2] * (*this)(l-1, m) - std::sqrt((l-1)*(l-1) - m*m) * r_squarenorm * (*this)(l-2, m)) / std::sqrt(l*l - m*m);
            elements.push_back(el.real());
            elements.push_back(el.imag());

            assert((*this)(l, m) == el);
            assert((*this)(l, -m) == (m%2 ? -std::conj(el) : std::conj(el)));
        }

        // case m = l-1
        el = r[2] * std::sqrt(2*l - 1) * (*this)(l-1, l-1);
        elements.push_back(el.real());
        elements.push_back(el.imag());
        assert((*this)(l, l - 1) == el);
        assert((*this)(l, -l+1) == (l%2 ? std::conj(el) : -std::conj(el)));

        // case m = l
        el = -sqrt((double)(2*l - 1)/(2*l)) * std::complex{r[0], r[1]} * (*this)(l-1, l-1);
        elements.push_back(el.real());
        elements.push_back(el.imag());       
        assert((*this)(l, l) == el);
        assert((*this)(l, -l) == (l%2 ? -std::conj(el) : std::conj(el)));
    }
}

std::complex<double> PowerSet<3>::operator()(unsigned l, int m) const
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
