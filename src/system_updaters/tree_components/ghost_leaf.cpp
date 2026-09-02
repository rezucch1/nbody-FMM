/**
 * @file ghost_leaf.cpp
 * @brief Implementation of GhostLeaf class methods.
 */

#include "ghost_leaf.hpp"

void GhostLeaf::calculateMC() {
    if (particles.empty()) return;

    mass_center = Tensor(dim);
    double total_mass = 0.0;

    for (const Particle *p : particles) {
        double m = p->get_weight();
        total_mass += m;
        mass_center += m * p->get_position();
    }

    if (total_mass > 0.0) {
        mass_center /= total_mass;
    }
}
