#include "power_set.hpp"

PowerSetI::PowerSetI(unsigned int L)
:BaseMultipoleSetI(L)
{ 
}

PowerSetI &PowerSetI::operator*=(double weight)
{
    for(auto &el : elements){
        el *= weight;
    }
    return *this;
}
