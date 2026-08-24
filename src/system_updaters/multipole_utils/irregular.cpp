#include "irregular.hpp"


Irregular::Irregular(unsigned int L, const Tensor &r)
: PowerSetI(L)
{
    elements.reserve((L+1)*(L+1));

    double r_squarenorm = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]; 

    elements[0] = 1/std::sqrt(r_squarenorm);
    elements[1]=r[2] / r_squarenorm;
    elements[2]=r[0] / r_squarenorm;
    elements[3]=r[1] / r_squarenorm;
    int i = 4;
    double mul_factor;

    for (int l=2; l<L+1; l++){
        for (int m=0; m<l-1; m++){
            mul_factor = 1 / (r_squarenorm * std::sqrt(l*l * m*m));
            elements[i] = mul_factor * (r[2]*(2*l-1)*elements[i - 2*l + 1] - elements[i - 4*(l+1)]);
            i++;
            elements[i] = mul_factor * (r[2]*(2*l-1)*elements[i - 2*l + 1] - elements[i - 4*(l+1)]);
            i++;
        }

        // case m = l-1
        mul_factor = std::sqrt(2*l - 1) / r_squarenorm;
        elements[i] = mul_factor * r[2]*elements[i - 2*l + 1];
        i++;
        elements[i] = mul_factor * r[2]*elements[i - 2*l + 1];
        i++;

        // case m = l
        mul_factor = std::sqrt(1 - 1/(2*l)) / r_squarenorm;
        elements[i] = mul_factor * (r[0]*elements[i-2*l-1] + r[1]*elements[i-2*l]);
        elements[i+1] = mul_factor * (r[0]*elements[i-2*l] - r[1]*elements[i-2*l-1]);
        i = i+2;        

    }

}

std::complex<double> Irregular::operator()(unsigned l, int m) const
{
    if (m = 0)
        return std::complex<double>(elements[l*l], 0);
    else if (m < 0)
        return std::complex<double>(elements[l*l + 2*m - 1], elements[l*l + 2*m]);
    else {
        if (m % 2 == 0)
            return std::complex<double>(elements[l*l + 2*m - 1], -elements[l*l + 2*m]);
        else
            return std::complex<double>(-elements[l*l + 2*m - 1], elements[l*l + 2*m]);
    }
}
