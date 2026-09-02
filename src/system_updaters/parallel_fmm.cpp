/**
 * @file parallel_fmm.cpp
 * @brief Implementation of ParallelFMM using ParallelTree.
 */

#include "parallel_fmm.hpp"
#include "tree_components/parallel_tree.hpp"

ParallelFMM::ParallelFMM() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
        MPI_Init(nullptr, nullptr);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
}

void ParallelFMM::update(const std::vector<Particle> &particles) {
    if (particles.empty()) return;

    // Encapsulate all parallel tree construction and parallel FMM calculation inside ParallelTree
    ParallelTree parallel_tree(particles.begin(), particles.end());
    parallel_tree.compute_parallel_accelerations();
}
