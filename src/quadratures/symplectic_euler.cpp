/**
 * @file symplectic_euler.cpp
 * @brief Implementation of Symplectic Euler integrator.
 */

#include "symplectic_euler.hpp"

void SyplecticEuler::integrate(SystemUpdateMethod *method, std::vector<Particle> &particles, const double &delta_t) {
  method->update(particles);
  for (unsigned int i = 0; i < particles.size(); ++i) {
    Particle &p = particles[i];
    velocity(p) += delta_t * p.get_acceleration();
    position(p) += delta_t * velocity(p);
  }
}