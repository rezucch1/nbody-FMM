/**
 * @file leaf.cpp
 * @brief Implementation of Leaf node P2M and acceleration evaluation methods.
 */

#include "leaf.hpp"
#include "multipole_set.hpp"

Leaf::Leaf(std::vector<std::vector<std::unique_ptr<NodeI>>> &allocator, unsigned int depth, unsigned int id_child, const Particle** particles_begin,
    const Particle** particles_end): NodeI(allocator, depth, id_child, particles_begin, particles_end){}

void Leaf::compute_multipoles(unsigned int L){
    calculateMC();
    multipole_set = std::make_unique<MultipoleSet<2>>(L);
    (*multipole_set) = 0;
    for( auto p = particles_begin; p< particles_end; ++p){
        Tensor r = (*p)->get_position() - mass_center;
        
        std::unique_ptr<PowerSetI> z;
        if (dim == 2)
            z = std::make_unique<PowerSet<2>>(L, r);
        else if (r.dim == 3){
            z = std::make_unique<PowerSet<3>>(L, r);
        }

        (*z) *= (*p)->get_weight();
        (*multipole_set) += *z;
    }
}

void Leaf::compute_acceleration(){
    for (auto i = particles_begin; i < particles_end; ++i){
        Tensor grad_i = - local_set->get_gradient((*i)->get_position() - mass_center);
        for (const auto &n : neighbours_list) if (n){
            for (auto j = ((Leaf*)n)->particles_begin; j < ((Leaf*)n)->particles_end; ++j){
                Tensor d = (*j)->get_position() - (*i)->get_position();
                grad_i -= (*j)->get_weight() * d / std::pow(d.squared_norm(), dim/2);
                if (dim % 2 == 1) grad_i /= d.norm();
            }
        }
    }
}

void Leaf::calculateMC()
{
    double total_mass = 0;

    // The mass_center should be initialized anyways, at least to 0;
    mass_center = 0 * (*particles_begin)->get_position();

    if (particles_begin == particles_end)
        return;

    for(auto p = particles_begin; p < particles_end; p++){
        total_mass += (*p)->get_weight();
        mass_center += (*p)->get_weight() * (*p)->get_position();
    }
    mass_center /= total_mass;

    return;
}
