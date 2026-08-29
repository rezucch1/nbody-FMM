/**
 * @file leaf.cpp
 * @brief Implementation of Leaf node P2M and acceleration evaluation methods.
 */

#include "leaf.hpp"
#include "multipole_set.hpp"

Leaf::Leaf(std::vector<std::unordered_map<uint64_t, std::unique_ptr<NodeI>>> &allocator, unsigned int depth, uint64_t id_child, std::vector<const Particle *>::const_iterator particles_begin, std::vector<const Particle *>::const_iterator particles_end)
    : NodeI(allocator, depth, id_child, (*particles_begin)->get_position().dim)
    , particles(particles_begin, particles_end){}

void Leaf::compute_multipoles(unsigned int L){
    calculateMC();
    if (dim == 2)
        multipole_set = std::make_unique<MultipoleSet<2>>(L);
    else if (dim == 3)
        multipole_set = std::make_unique<MultipoleSet<3>>(L);

    (*multipole_set) = 0;
    for(const auto &p : particles){
        Tensor r = p->get_position() - mass_center;
        
        std::unique_ptr<PowerSetI> z;
        if (dim == 2)
            z = std::make_unique<PowerSet<2>>(L, r);
        else if (dim == 3)
            z = std::make_unique<PowerSet<3>>(L, r);

        (*z) *= p->get_weight();
        (*multipole_set) += *z;
    }
    if ((*multipole_set)(0).real() == 0)
        printf("break\n");
    for(int n = 0; n <= L; ++n)
        for (int m = 0; m <= n; ++m)
            if (std::isnan((*multipole_set)(n, m).real()) || std::isnan((*multipole_set)(n, m).imag()))
                printf("break\n");
}

void Leaf::compute_acceleration(){
    for (const auto &i : particles){
        Tensor grad_i(dim);
        grad_i *= 0;
        if (local_set) {
            grad_i = local_set->get_gradient(i->get_position() - mass_center);
        }

        // 1. Direct P2P interactions with other particles in the SAME leaf cell
        for (const auto &j : particles) {
            if (j != i) {
                Tensor d = j->get_position() - i->get_position();
                grad_i -= j->get_weight() * d / std::pow(d.squared_norm(), (dim) / 2);
            }
        }

        // 2. Direct P2P interactions with particles in NEIGHBOR leaf cells
        for (const auto &n : neighbours_list) if (n){
            for (const auto &j : ((Leaf*)n)->particles){
                Tensor d = j->get_position() - i->get_position();
                grad_i -= j->get_weight() * d / std::pow(d.squared_norm(), (dim) / 2);
            }
        }
        i->compute_new_accelaration(grad_i);
    }
}

void Leaf::calculateMC()
{
    double total_mass = 0;

    // The mass_center should be initialized anyways, at least to 0;
    mass_center = 0 * (*particles.begin())->get_position();

    if (particles.begin() == particles.end())
        return;

    for(const auto p : particles){
        total_mass += p->get_weight();
        mass_center += p->get_weight() * p->get_position();
    }
    mass_center /= total_mass;

    return;
}
