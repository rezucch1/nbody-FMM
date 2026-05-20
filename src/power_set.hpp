#pragma once
#include <vector>
#include "tensor.hpp"
#include <complex>

class PowerSetI{
    protected :
    std::vector<double> elements;
    unsigned int L;
    PowerSetI(unsigned int L);
    public :
    virtual std::complex<double> operator()(unsigned l, int m)=0;
    PowerSetI &operator *= (double mass);

    friend MultipoleSetI &MultipoleSetI::operator+=(const PowerSetI &);
};

#include "multipole_set.hpp"

template<unsigned int dim>
class PowerSet : public PowerSetI{

};

#include "power_set2D.hpp"
#include "power_set3D.hpp" 


