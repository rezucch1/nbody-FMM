#pragma once
#include "base_multipole_set_i.hpp"
#include "math_utils/tensor.hpp"

class LocalSetI : public BaseMultipoleSetI{
  public:
    virtual LocalSetI &operator+=(const LocalSetI *other) = 0;
    virtual LocalSetI* distribute_parent_with_distance (const Tensor &d) const=0;
  protected:
    LocalSetI(unsigned int L) : BaseMultipoleSetI(L) {};
};

template <unsigned int dim>
class LocalSet : public LocalSetI{

};

#include "local_set2D.hpp"
#include "local_set3D.hpp"