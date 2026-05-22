#include "multipole_set.hpp"

MultipoleSetI::MultipoleSetI(unsigned int L)
:BaseMultipoleSetI(L)
{ 

}

MultipoleSetI &MultipoleSetI::operator+=(const PowerSetI &z)
{
    for(int i=0; i<elements.size(); ++i){
        this->elements[i] += z.elements[i];
    }
}

MultipoleSetI &MultipoleSetI::operator+=(const MultipoleSetI &z)
{
    
    for(int i=0; i<elements.size(); ++i){
        this->elements[i] += z.elements[i];
    }
}

const std::vector<double> &MultipoleSetI::getElements(const MultipoleSetI &_this)
{
    return _this.elements;
}
