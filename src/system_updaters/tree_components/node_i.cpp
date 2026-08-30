/**
 * @file node_i.cpp
 * @brief Implementation of NodeI base methods, Morton interaction arithmetic, M2L collection, and L2L propagation.
 */

#include "node_i.hpp"
#include "node.hpp"

Node &NodeI::get_parent(){
  return *(Node*)allocator[depth - 1][id_child >> dim].get();
}

const std::unique_ptr<MultipoleSetI> &NodeI::getMultipoleSet(const NodeI &_this)
{
    return _this.multipole_set;
}

const Tensor &NodeI::getMassCenter(const NodeI &_this)
{
    return _this.mass_center;
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

            auto it = allocator[depth].find(neighbour_id);
            if (it != allocator[depth].end() && it->second) {
                NodeI* neighbour_node = it->second.get();
                if (is_neighbour) {
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
    for (; s < interaction_list.cend() && (*s == nullptr || !getMultipoleSet(**s)); ++s);
    if (s == interaction_list.cend()){
        local_set =  nullptr;
        return;
    }
    local_set = getMultipoleSet(**s)->to_local(mass_center - getMassCenter(**s));

    ++s;

    for (; s != interaction_list.cend(); ++s) if (*s && getMultipoleSet(**s)){
        *local_set += getMultipoleSet(**s)->to_local(mass_center - getMassCenter(**s)).get();
    }
}

void NodeI::propagate_locals(const LocalSetI *parent_local){
    if (parent_local){
        if (local_set){
            *local_set += parent_local;
        }
        else
            local_set = parent_local->clone();
    }
}
