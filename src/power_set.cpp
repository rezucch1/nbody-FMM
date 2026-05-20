#include "power_set.hpp"

PowerSetI::PowerSetI(unsigned int L)
:L(L)
{ 

}

PowerSetI &PowerSetI::operator*=(double mass)
{
    for(auto &el : elements){
        el *= mass;
    }
    return *this;
}
