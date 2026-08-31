#pragma once

#include "system_update_method.hpp"
#include "tree_components/tree.hpp"
#include <mpi.h>
#include <vector>
#include <cstdint>

/**
 * @file parallel_fmm.hpp
 * @brief Header for MPI Parallel Fast Multipole Method (ParallelFMM) solver.
 */

/**
 * @class ParallelFMM
 * @brief MPI-parallelized FMM solver using coarse branch Morton domain decomposition.
 */
class ParallelFMM : public SystemUpdateMethod {
  public:
    /**
     * @brief Constructs ParallelFMM solver and initializes MPI rank parameters.
     */
    ParallelFMM();

    /**
     * @brief Virtual destructor.
     */
    virtual ~ParallelFMM() = default;

    /**
     * @brief Computes acceleration vectors in parallel using MPI across ranks.
     * @param particles Reference to global system particles array.
     */
    virtual void update(const std::vector<Particle> &particles) override;

    /**
     * @brief Gets current MPI process rank.
     * @return Rank index.
     */
    int get_rank() const { return rank; }

    /**
     * @brief Gets total number of MPI processes.
     * @return Number of ranks.
     */
    int get_num_ranks() const { return num_ranks; }

  private:
    int rank;
    int num_ranks;
};
