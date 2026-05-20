#pragma once
#include <vector>
#include "tensor.hpp"
#include <complex>

class MultipoleSetI;
#include "power_set.hpp"


class MultipoleSetI{
    protected :
    std::vector<double> elements;
    unsigned int L;
MultipoleSetI(unsigned int L);
    public :
    virtual std::complex<double> operator()(unsigned l, int m)=0;
    MultipoleSetI &operator += (const PowerSetI &z);
};

template<unsigned int dim>
class MultipoleSet : public   MultipoleSetI{

};

#include "multipole_set2D.hpp"
#include "multipole_set3D.hpp" 


