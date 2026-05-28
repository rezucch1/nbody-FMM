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
    return *this;
}

MultipoleSetI &MultipoleSetI::operator+=(const MultipoleSetI &z)
{
    
    for(int i=0; i<elements.size(); ++i){
        this->elements[i] += z.elements[i];
    }
    return *this;
}

MultipoleSetI &MultipoleSetI::operator=(std::nullptr_t)
{
    int i = 0;
    for (; i < elements.size(); ++i)
        elements[i] = 0.0;
    for (; elements.capacity(); ++i)
        elements.push_back(0.0);
    return *this;
}

const std::vector<double> &MultipoleSetI::getElements(const MultipoleSetI &_this)
{
    return _this.elements;
}
