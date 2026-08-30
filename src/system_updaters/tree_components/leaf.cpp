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
void Leaf::compute_acceleration(){
    for (const auto &i : particles){

        Tensor grad_i;

        // ============================================================
        // 2D
        // ============================================================
        if (dim == 2) {

            grad_i = Tensor{0.0, 0.0};

            // FMM: contributo del local expansion
            if (local_set) {
                Tensor d = i->get_position() - mass_center;

                grad_i += local_set->get_gradient(d);

            }

            // P2P: particelle nella stessa leaf
            for (const auto &j : particles) {
                if (j != i) {
                    Tensor d = i->get_position() - j->get_position();

                    // Kernel: -log(r)
                    grad_i -= j->get_weight()
                            * d / d.squared_norm();
                }
            }

            // P2P: particelle nelle leaf vicine
            for (const auto &n : neighbours_list) {
                if (n) {
                    for (const auto &j : ((Leaf*)n)->particles) {
                        if (i != j) {
                            Tensor d = i->get_position() - j->get_position();

                            // Kernel: -log(r)
                            grad_i -= j->get_weight()
                                    * d / d.squared_norm();
                        }
                    }
                }
            }
        }

        // ============================================================
        // 3D
        // ============================================================
        else if (dim == 3) {

            grad_i = Tensor{0.0, 0.0, 0.0};

            // FMM: contributo del local expansion
            if (local_set) {
                Tensor d = i->get_position() - mass_center;

                grad_i += local_set->get_gradient(d);

            }

            // P2P: particelle nella stessa leaf
            for (const auto &j : particles) {
                if (j != i) {
                    Tensor d = i->get_position() - j->get_position();

                    // Kernel: 1/r
                    grad_i -= j->get_weight()
                            * d / std::pow(d.squared_norm(), 1.5);
                }
            }

            // P2P: particelle nelle leaf vicine
            for (const auto &n : neighbours_list) {
                if (n) {
                    for (const auto &j : ((Leaf*)n)->particles) {
                        if (i != j) {
                            Tensor d = i->get_position() - j->get_position();

                            // Kernel: 1/r
                            grad_i -= j->get_weight()
                                    * d / std::pow(d.squared_norm(), 1.5);
                        }
                    }
                }
            }
        }

        // Aggiorna l'accelerazione
        i->compute_new_accelaration(grad_i);
    }
}
