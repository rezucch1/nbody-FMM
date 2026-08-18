#include "node_i.hpp"

const std::unique_ptr<MultipoleSetI> &NodeI::getMultipoleSet(const NodeI &_this)
{
    return _this.multipole_set;
}

const Tensor &NodeI::getMassCenter(const NodeI &_this)
{
    return _this.mass_center;
}

void NodeI::get_partition(std::vector<std::tuple<Particle *, int, int>> &partitions, int level, int partition_id) const{
    std::vector<std::tuple<Particle *, int, int>> sub_partition;
    sub_partition.reserve(particles_end - particles_begin);

    for (Particle **p = particles_begin; p < particles_end; ++p)
        sub_partition.push_back(std::make_tuple(*p, level, partition_id));
    
    partitions.insert(partitions.end(), sub_partition.begin(), sub_partition.end());
}
