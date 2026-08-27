/**
 * @file runge_kutta_4.cpp
 * @brief Implementation of classical 4th-order Runge-Kutta integrator.
 */

#include "runge_kutta_4.hpp"
#include <cassert>
#include <cmath>

void RungeKutta4::integrate(SystemUpdateMethod *method, std::vector<Particle> &particles, const double &delta_t){
  // a_k1 + 2a_k2 + 2a_k3 + a_k4
  std::unique_ptr<Tensor[]> accelerations_sum = std::make_unique<Tensor[]>(particles.size());

  // v_k1 + 2v_k2 + 2v_k3 + v_k4
  std::unique_ptr<Tensor[]> velocity_sum = std::make_unique<Tensor[]>(particles.size());

  std::vector<Particle> particles_k;
  particles_k.reserve(particles.size());

  std::unique_ptr<Tensor[]> velocity_k = std::make_unique<Tensor[]>(particles.size());

  // stage 1
  std::copy(particles.begin(), particles.end(), std::back_inserter(particles_k));

  for (unsigned int i = 0; i < particles.size(); ++i){
    velocity_k[i] = particles[i].get_velocity();
  }

  method->update(particles_k);
  for (unsigned int i = 0; i < particles.size(); ++i){
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
