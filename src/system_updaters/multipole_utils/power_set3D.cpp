/**
 * @file power_set3D.cpp
 * @brief Implementation of 3D PowerSet regular solid harmonics evaluation.
 */

#include "power_set3D.hpp"
#include <cmath>

PowerSet<3>::PowerSet(unsigned int L, const Tensor &r)
: PowerSetI(L)
{
    elements.reserve((L+1)*(L+1));
    elements.push_back(1);
    elements.push_back(r[2]); // z
    elements.push_back(-1/M_SQRT2 * r[0]);
    elements.push_back(-1/M_SQRT2 * r[1]);

    double r_squarenorm = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
    int i = 4;
    double mul_factor;

    for (int l=2; l<=L; l++){
        elements.push_back(mul_factor * (r[2]*(2*l-1)*elements[i - 2*l + 1] - (l - 1) * r_squarenorm * elements[i - 4*(l-1)]));
        ++i;
        for (int m=1; m<l-1; m++){
            mul_factor = 1 / std::sqrt(l*l - m*m);
            elements.push_back(mul_factor * (r[2]*(2*l-1)*elements[i - 2*l + 1] - std::sqrt((l + m -1) * (l - m - 1)) * r_squarenorm * elements[i - 4*(l-1)]));
            i++;
            elements.push_back(mul_factor * (r[2]*(2*l-1)*elements[i - 2*l + 1] - std::sqrt((l + m -1) * (l - m - 1)) * r_squarenorm * elements[i - 4*(l-1)]));
            i++;
        }

        // case m = l-1
        mul_factor = std::sqrt(2*l - 1);
        elements.push_back(mul_factor * r[2]*elements[i - 2*l + 1]);
        i++;
        elements.push_back(mul_factor * r[2]*elements[i - 2*l + 1]);
        i++;

        // case m = l
        mul_factor = std::sqrt(1 - 1/(2*l));
        elements.push_back(mul_factor * (r[0]*elements[i-2*l-1] + r[1]*elements[i-2*l]));
        elements.push_back(mul_factor * (r[0]*elements[i-2*l] - r[1]*elements[i-2*l-1]));
        i = i+2;        
    }
}

std::complex<double> PowerSet<3>::operator()(unsigned l, int m) const
{
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
