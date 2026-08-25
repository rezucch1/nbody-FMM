#include "node_i.hpp"
#include "node.hpp"

Node &NodeI::get_parent(){
  return *(Node*)allocator[depth - 1][id_child / 4].get();
}

const std::unique_ptr<MultipoleSetI> &NodeI::getMultipoleSet(const NodeI &_this)
{
    return _this.multipole_set;
}

const Tensor &NodeI::getMassCenter(const NodeI &_this)
{
    return _this.mass_center;
}

void NodeI::get_partition(std::vector<std::tuple<Particle *, int, int>> &partitions) const{
    std::vector<std::tuple<Particle *, int, int>> sub_partition;
    sub_partition.reserve(particles_end - particles_begin);

    for (Particle **p = particles_begin; p < particles_end; ++p)
        sub_partition.push_back(std::make_tuple(*p, depth, id_child));
    
    partitions.insert(partitions.end(), sub_partition.begin(), sub_partition.end());
}

void NodeI::compute_interaction_list(){
    neighbours_list.reserve(std::pow(3, dim) - 1);
    interaction_list.reserve(std::pow(3, dim) * ((1<<dim) -1));
    
    int cartesian_coor[dim];
    for (unsigned int i = 0; i < dim; ++i)
        cartesian_coor[i] = 0;
    
    unsigned int id_child_temp = id_child;
    unsigned int digit = 0;
    while (id_child_temp > 0){
        for (unsigned int i = 0; i < dim; ++i){
            cartesian_coor[i] += (id_child_temp % 2) << digit;
            id_child_temp /= 2;
        }
        ++digit;
    }

    int interaction_start[dim];
    int interaction_end[dim];
    int direction[dim];

    for (int i = 0; i < dim; ++i){
        interaction_start[i] = - (cartesian_coor[i] % 2) - 2;
        interaction_end[i] = 3 - cartesian_coor[i] % 2;

        direction[i] = interaction_start[i];
    }

    bool is_loop_ended;

    do{
        int n_cartesian_coor[dim];
        bool is_out = false;
        for(int i = 0; i < dim && !is_out; i++){
            n_cartesian_coor[i] = cartesian_coor[i] + direction[i];
            if (n_cartesian_coor[i] < 0) is_out = true;
        }

        if (!is_out){
            unsigned int neighbour_id = 0;
            for (unsigned int i = 0; i < dim; ++i){
                digit = i;
                while (n_cartesian_coor[i] > 0){
                    neighbour_id += (n_cartesian_coor[i] % 2) << digit;
                    n_cartesian_coor[i] /= 2;
                    digit += dim;
                }
            }

            bool is_neighbour = true;
            for(unsigned int i = 0; i < dim && is_neighbour; ++i)
                if (direction[i] < -1 || direction[i] > 1)
                    is_neighbour = false;

            NodeI* neighbour_node;
            if (neighbour_id >= 0 && neighbour_id < allocator[depth].size()){
                neighbour_node = allocator[depth][neighbour_id].get();
                // else 
                //     neighbour_node = nullptr;

                if (is_neighbour){
                    if (neighbour_id != id_child)
                        neighbours_list.push_back(neighbour_node);
                } else
                    interaction_list.push_back(neighbour_node);
            }
        }

        is_loop_ended = true;
        for(int i = 0; i < dim && is_loop_ended; ++i){
            if (direction[i] == interaction_end[i]){
                direction[i] = interaction_start[i];
            } else {
                ++direction[i];
                is_loop_ended = false;
            }
        }
    }while (!is_loop_ended);
    
}

unsigned int NodeI::get_id() const{
  return id_child;
}

void NodeI::collect_multipoles_to_locals(){
    auto s = interaction_list.cbegin();
    for (; s < interaction_list.cend(); ++s);
    if (s == interaction_list.cend()){
        local_set = nullptr;
        return;
    }
    local_set = getMultipoleSet(**s)->to_local(mass_center - getMassCenter(**s));

    ++s;

    for (; s < interaction_list.cend(); ++s){
        *local_set += getMultipoleSet(**s)->to_local(mass_center - getMassCenter(**s)).get();
    }
}

void NodeI::propagate_locals(LocalSetI *parent_local){
    if (parent_local){
        if (local_set){
            *local_set += parent_local;
            delete parent_local;
        }
        else
            local_set = std::unique_ptr<LocalSetI>(parent_local);
    }
}
