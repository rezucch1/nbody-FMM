#include "naive_update.hpp"
#include <cmath>

const Tensor* NaiveUpdate::update(const std::vector<Particle> &particles) {
  Tensor particle_acceleration;
  next_acceleration.clear();
  for (const Particle &p : particles){
    particle_acceleration = 0 * p.get_position();
    for (const Particle &other : particles) if (&p != &other){
      Tensor distance = p.get_position() - other.get_position();
      double squared_distance_norm = distance.squared_norm();
      particle_acceleration += other.get_weight() / ((squared_distance_norm + 1.0e-4) * std::sqrt(squared_distance_norm)) * distance ;
    }
    particle_acceleration *= -GRAV_CONST;
    next_acceleration.push_back(particle_acceleration);
  }
  return next_acceleration.data();
}