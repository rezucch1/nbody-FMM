#include "multipole_set3D.hpp"


inline MultipoleSet<3>::MultipoleSet(unsigned int L)
: MultipoleSetI(L)
{
}

std::complex<double> MultipoleSet<3>::operator()(unsigned l, int m) const
{
    return std::complex<double>();
}
