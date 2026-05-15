#pragma once
#include <vector>
#include "tensor.hpp"


class MultipoleSetI{
    protected :
    std::vector<double> elements;
    unsigned int L;
    MultipoleSetI(unsigned int L);
};

template<unsigned int dim>
class MultipoleSet : public MultipoleSetI{

};
