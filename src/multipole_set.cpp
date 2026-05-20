#include "multipole_set.hpp"

MultipoleSetI::MultipoleSetI(unsigned int L)
:L(L)
{ 

}

MultipoleSetI &MultipoleSetI::operator+=(const PowerSetI &z)
{
    for(int i=0; i<elements.size(); ++i){
        this->elements[i] += z.elements[i];
    }
}
