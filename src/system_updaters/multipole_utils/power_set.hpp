#pragma once
#include <vector>
#include "math_utils/tensor.hpp"
#include <complex>
#include "base_multipole_set_i.hpp"

class PowerSetI;
#include "multipole_set.hpp"

class PowerSetI : public BaseMultipoleSetI{
    protected :
    PowerSetI(unsigned int L);
    public :
    PowerSetI &operator *= (double mass);
};

template<unsigned int dim>
class PowerSet : public PowerSetI{

};

#include "power_set2D.hpp"
// #include "power_set3D.hpp" 


