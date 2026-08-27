/**
 * @file multipole_set.cpp
 * @brief Implementation of MultipoleSetI base operators.
 */

#include "multipole_set.hpp"

MultipoleSetI::MultipoleSetI(unsigned int L)
:BaseMultipoleSetI(L)
{ 

}

MultipoleSetI &MultipoleSetI::operator+=(const PowerSetI &z)
{
    for(size_t i=0; i<elements.size(); ++i){
        this->elements[i] += z.elements[i];
    }
    return *this;
}

MultipoleSetI &MultipoleSetI::operator+=(const MultipoleSetI &z)
{
    for(size_t i=0; i<elements.size(); ++i){
        this->elements[i] += z.elements[i];
    }
    return *this;
}

MultipoleSetI &MultipoleSetI::operator=(std::nullptr_t)
{
    size_t i = 0;
    for (; i < elements.size(); ++i)
        elements[i] = 0.0;
    for (; i < elements.capacity(); ++i)
        elements.push_back(0.0);
    return *this;
}

const std::vector<double> &MultipoleSetI::getElements(const MultipoleSetI &_this)
{
    return _this.elements;
}

std::ostream &operator<< (std::ostream &os, const MultipoleSetI &multipole){
    return os;
}
