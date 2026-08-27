/**
 * @file FMM.cpp
 * @brief Implementation of Fast Multipole Method (FMM) update solver.
 */

#include "FMM.hpp"
#include "tree_components/tree.hpp"

void FMM::update(const std::vector<Particle> &particles) {
  Tree tree(particles.begin(), particles.end());
  return tree.compute_accelerations();
}