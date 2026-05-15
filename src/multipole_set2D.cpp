#include "multipole_set2D.hpp"


inline MultipoleSet<2>::MultipoleSet(unsigned int L, const Tensor &r)
: MultipoleSetI(L)
{
    elements.reserve(2*L+1);
    elements[0] = 1; 
    
}

