#include "power_set3D.hpp"


inline PowerSet<3>::PowerSet(unsigned int L, const Tensor &r)
: PowerSetI(L)
{
    elements.reserve((L+1)*(L+1));
    elements[0] = 1;
    elements[1]=r[2];
    elements[2]=r[0];
    elements[3]=r[1];

    double r_squarenorm = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]; 
    int i = 4;

    for (int l=2; l<L+1; l++){
        for (int m=0; m<l-1; m++){
            elements[i] = 1/((l*l)-(m*m)) * (r[2]*(2*l-1)*elements[i - 2*l + 1] - r_squarenorm * elements[i - 4*(l+1)]);
            i++;
            elements[i] = 1/((l*l)-(m*m)) * (r[2]*(2*l-1)*elements[i - 2*l + 1] - r_squarenorm * elements[i - 4*(l+1)]);
            i++;
        }

        // case m = l-1
        elements[i] = 1/(2*l - 1) * r[2]*(2*l-1)*elements[i - 2*l + 1];
        i++;
        elements[i] = 1/(2*l - 1) * r[2]*(2*l-1)*elements[i - 2*l + 1];
        i++;

        // case m = l
        elements[i] = 1/(2*l) * (r[0]*elements[i-2*l-1] - r[1]*elements[i-2*l]);
        elements[i+1] = 1/(2*l) * (r[0]*elements[i-2*l] + r[1]*elements[i-2*l-1]);
        i = i+2;        

    }

}

std::complex<double> PowerSet<3>::operator()(unsigned l, int m)
{
    return std::complex<double>();
}
