#include "symplectic_euler.hpp"

void SyplecticEuler::integrate(SystemUpdateMethod *method, std::vector<Particle> &particles, const double &delta_t){
  const Tensor* acceleration = method->update(particles);
  for (unsigned int i = 0; i < particles.size(); ++i){
    Particle &p = particles[i];
    velocity(p) += delta_t * acceleration[i];
    position(p) += delta_t * velocity(p);
  }
}