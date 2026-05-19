#pragma once
#include <vector>
#include "tensor.hpp"
#include <complex>


class MultipoleSetI{
    protected :
    std::vector<double> elements;
    unsigned int L;
    MultipoleSetI(unsigned int L);
    public :
    virtual std::complex<double> operator()(unsigned l, int m)=0;
};

template<unsigned int dim>
class MultipoleSet : public MultipoleSetI{

};


