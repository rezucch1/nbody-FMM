#pragma once
#include <vector>
#include "math_utils/tensor.hpp"
#include <complex>
#include "base_multipole_set_i.hpp"
#include <iostream>

class MultipoleSetI;
#include "power_set.hpp"

class MultipoleSetI : public BaseMultipoleSetI{
    protected :
    MultipoleSetI(unsigned int L);

    virtual std::ostream &output_to_stream (std::ostream &os) const = 0;
    friend std::ostream &operator<<(std::ostream &os, const MultipoleSetI &multipole);

    public :
    MultipoleSetI &operator += (const PowerSetI &z);
    MultipoleSetI &operator += (const MultipoleSetI &z);
    MultipoleSetI &operator = (std::nullptr_t);
    virtual std::unique_ptr<MultipoleSetI> weigh_children_with_distance (const Tensor &d) const=0;
    static const std::vector<double> &getElements(const MultipoleSetI &_this);
};

std::ostream &operator<< (std::ostream &os, const MultipoleSetI &multipole);

template<unsigned int dim>
class MultipoleSet : public   MultipoleSetI{

};

#include "multipole_set2D.hpp"
// #include "multipole_set3D.hpp" 


