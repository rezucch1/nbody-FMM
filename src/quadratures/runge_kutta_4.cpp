/**
 * @file runge_kutta_4.cpp
 * @brief Implementation of classical 4th-order Runge-Kutta integrator.
 */

#include "runge_kutta_4.hpp"
#include <cassert>
#include <cmath>

void RungeKutta4::integrate(SystemUpdateMethod *method, std::vector<Particle> &particles, const double &delta_t){
  // a_k1 + 2a_k2 + 2a_k3 + a_k4
  std::vector<Tensor> accelerations_sum;
  // REVIEW: these vectors are indexed from 0 below after reserve(), but reserve()
  // does not create elements. Construct/resize all three vectors to particles.size()
  // before indexing; the present code has undefined behaviour on every non-empty step.
  accelerations_sum.reserve(particles.size());

  // v_k1 + 2v_k2 + 2v_k3 + v_k4
  std::vector<Tensor> velocity_sum;
  velocity_sum.reserve(particles.size());

  std::vector<Particle> particles_k;
  particles_k.reserve(particles.size());

  std::vector<Tensor> velocity_k;
  velocity_k.reserve(particles.size());

  // stage 1
  std::copy(particles.begin(), particles.end(), std::back_inserter(particles_k));

  for (unsigned int i = 0; i < particles.size(); ++i){
    velocity_k[i] = particles[i].get_velocity();
  }

  method->update(particles_k);
  for (unsigned int i = 0; i < particles.size(); ++i){
    // acceleration_sum has zero size. you have to use emplace_back!
    accelerations_sum[i] = particles_k[i].get_acceleration();
    assert(!std::isnan(accelerations_sum[i].squared_norm()));
    velocity_sum[i] = velocity_k[i];
    assert(&velocity_sum[i] != &velocity(particles[i]));
  }

  // stage 2
  std::copy(particles.begin(), particles.end(), particles_k.begin());
  for (unsigned int i = 0; i < particles.size(); ++i){
    position(particles_k[i]) += delta_t / 2 * velocity_k[i];
  }

  for (unsigned int i = 0; i < particles.size(); ++i){
    velocity_k[i] = particles[i].get_velocity() + delta_t / 2 * particles_k[i].get_acceleration();
  }

  method->update(particles_k);
  for (unsigned int i = 0; i < particles.size(); ++i){
    accelerations_sum[i] += 2*particles_k[i].get_acceleration();
    velocity_sum[i] += 2*velocity_k[i];
  }

  // stage 3
  std::copy(particles.begin(), particles.end(), particles_k.begin());
  for (unsigned int i = 0; i < particles.size(); ++i){
    position(particles_k[i]) += delta_t / 2 * velocity_k[i];
  }

  for (unsigned int i = 0; i < particles.size(); ++i){
    velocity_k[i] = particles[i].get_velocity() + delta_t / 2 * particles_k[i].get_acceleration();
  }

  method->update(particles_k);
  for (unsigned int i = 0; i < particles.size(); ++i){
    accelerations_sum[i] += 2*particles_k[i].get_acceleration();
    velocity_sum[i] += 2*velocity_k[i];
  }

  // stage 4
  std::copy(particles.begin(), particles.end(), particles_k.begin());
  for (unsigned int i = 0; i < particles.size(); ++i){
    position(particles_k[i]) += delta_t * velocity_k[i];
  }
  for (unsigned int i = 0; i < particles.size(); ++i){
    velocity_k[i] = particles[i].get_velocity() + delta_t * particles_k[i].get_acceleration();
  }

  method->update(particles_k);

  for (unsigned int i = 0; i < particles.size(); ++i){
    accelerations_sum[i] += particles_k[i].get_acceleration();
    velocity_sum[i] += velocity_k[i];
  }

  // final next velocity and position
  for (unsigned int i = 0; i < particles.size(); ++i){
    velocity(particles[i]) += delta_t / 6 * accelerations_sum[i];
    position(particles[i]) += delta_t / 6 * velocity_sum[i];
  }
}
