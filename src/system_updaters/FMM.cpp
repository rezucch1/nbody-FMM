/**
 * @file FMM.cpp
 * @brief Implementation of Fast Multipole Method (FMM) update solver.
 */

#include "FMM.hpp"
#include "tree_components/tree.hpp"

const Tensor* FMM::update(const std::vector<Particle> &particles) {
  next_acceleration.clear();
  Tree tree(particles.begin(), particles.end());
  return tree.get_accelerations(next_acceleration, particles.data(), particles.size());
}