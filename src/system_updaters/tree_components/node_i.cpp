#include "node_i.hpp"

const std::unique_ptr<MultipoleSetI> &NodeI::getMultipoleSet(const NodeI &_this)
{
    return _this.multipole_set;
}

const Tensor &NodeI::getMassCenter(const NodeI &_this)
{
    return _this.mass_center;
}
