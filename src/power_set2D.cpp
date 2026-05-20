#include "power_set2D.hpp"
#include "math_utils.hpp"

inline PowerSet<2>::PowerSet(unsigned int L, const Tensor &r)
: PowerSetI(L)
{
    elements.reserve(2*L+1);
    elements[0] = 1;
    elements[1] = r[0];
    elements[2] = r[1];
    for(int i=3; i<2*L+1; i=i+2){
        std::tie(elements[i],elements[i+1])=mult_complex(elements[i-2], elements[i-1], r[0], r[1]);
        //elements[i] = elements[i-2]*r[0] - elements[i-1]*r[1];
        //elements[i+1] = elements[i-2]*r[1] + elements[i-1]*r[0];
    } 
}

std::complex<double> PowerSet<2>::operator()(unsigned l, int m)
{
    if (l=0) return std::complex<double>(elements[0], 0);
    return std::complex<double>(elements[2*l-1], elements[2*l]);
}
