#include "multipole_set2D.hpp"


inline MultipoleSet<2>::MultipoleSet(unsigned int L)
: MultipoleSetI(L)
{
}

std::complex<double> MultipoleSet<2>::operator()(unsigned l, int m)
{
    if (l=0) return std::complex<double>(elements[0], 0);
    return std::complex<double>(elements[2*l-1], elements[2*l]);
}
