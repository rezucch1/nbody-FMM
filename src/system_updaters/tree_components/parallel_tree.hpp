#pragma once

#include "tree.hpp"
#include <mpi.h>
#include <vector>
#include <cstdint>
#include <unordered_map>

/**
 * @file parallel_tree.hpp
 * @brief Header for ParallelTree class performing MPI parallel tree partitioning and FMM passes.
 */

/**
 * @class ParallelTree
 * @brief Inherits from Tree to perform MPI-parallel domain bounds reduction, branch ownership filtering, and parallel FMM calculation.
 */
class ParallelTree : public Tree {
  public:
    /**
     * @brief Constructs a ParallelTree for a particle range across MPI ranks.
     * @tparam InputIt Iterator over Particle pointers/objects.
     * @param begin Start iterator.
     * @param end End iterator.
     */
    template <typename InputIt>
    ParallelTree(InputIt begin, InputIt end) : Tree() {
        int initialized = 0;
        MPI_Initialized(&initialized);
        if (!initialized) {
            rank = 0;
            num_ranks = 1;
        } else {
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
        }

        init_parallel_tree(&*begin, &*end);
    }

    /**
     * @brief Destructor.
     */
    virtual ~ParallelTree() = default;

    /**
     * @brief Executes parallel FMM acceleration pass.
     */
    void compute_parallel_accelerations();

    /**
     * @brief Scatters coarse branch multipoles at d_coarse across all MPI ranks using MPI_Allgather.
     * @param L Maximum degree limit L of multipole expansions.
     */
    void scatter_coarse_multipoles(unsigned int L);

    /**
     * @brief Scatters deep ghost node multipoles (d > d_coarse) near process boundaries across ranks.
     * @param L Maximum degree limit L of multipole expansions.
     */
    void scatter_deep_ghost_multipoles(unsigned int L);

    int get_rank() const { return rank; }
    int get_num_ranks() const { return num_ranks; }

  protected:
    /**
     * @brief Initializes parallel tree hierarchy using MPI bounds reduction and coarse branch filtering.
     * @param begin Start particle pointer.
     * @param end End particle pointer.
     */
    void init_parallel_tree(const Particle *begin, const Particle *end);

  private:
    int rank;
    int num_ranks;
    unsigned int dim;
    unsigned int d_coarse;
};
