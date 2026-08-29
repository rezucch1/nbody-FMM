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
    if (dim == 2)
        multipole_set = std::make_unique<MultipoleSet<2>>(L);
    else if (dim == 3)
        multipole_set = std::make_unique<MultipoleSet<3>>(L);

    (*multipole_set) = 0;
    for( auto p = particles_begin; p< particles_end; ++p){
        Tensor r = (*p)->get_position() - mass_center;
        
        std::unique_ptr<PowerSetI> z;
        if (dim == 2)
            z = std::make_unique<PowerSet<2>>(L, r);
        else if (dim == 3)
            z = std::make_unique<PowerSet<3>>(L, r);

        (*z) *= (*p)->get_weight();
        (*multipole_set) += *z;
    }
    if ((*multipole_set)(0).real() == 0)
        printf("break\n");
}

void Leaf::compute_acceleration(){
    for (auto i = particles_begin; i < particles_end; ++i){
        Tensor grad_i = local_set->get_gradient((*i)->get_position() - mass_center);
        
        for (const auto &n : neighbours_list) if (n){
            for (auto j = ((Leaf*)n)->particles_begin; j < ((Leaf*)n)->particles_end; ++j){
                if (*i == *j)
                    continue;
                
                Tensor d = (*i)->get_position() - (*j)->get_position();

                 if (dim == 2) {
                    // Kernel: -log(r)
                    grad_i -= (*j)->get_weight()
                            * d / d.squared_norm();
                }
                else if (dim == 3) {
                    // Kernel: 1/r
                    grad_i -= (*j)->get_weight()
                            * d / std::pow(d.squared_norm(), 1.5);
                }

            }
        }
        (*i)->compute_new_accelaration(grad_i);
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
