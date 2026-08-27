#pragma once
#include <memory>
#include <vector>
#include <tuple>

#include "node.hpp"

/**
 * @file tree.hpp
 * @brief Tree class managing spatial decomposition, Morton ordering, and FMM tree execution.
 */

/**
 * @class Tree
 * @brief Top-level tree data structure performing iterative spatial partitioning and FMM execution.
 */
class Tree {
  public:
    /**
     * @brief Constructs FMM tree over a particle iterator range.
     * @tparam InputIt Iterator type pointing to Particle.
     * @param begin Start iterator.
     * @param end End iterator.
     */
    template <class InputIt>
    Tree(InputIt begin, InputIt end) { init_tree(&*begin, &*end); };

    /**
     * @brief Computes gravitational acceleration vectors for all particles using FMM.
     * @param acceleration_vector Reference to vector storing computed accelerations.
     * @param begin Pointer to start of original particle array.
     * @param size Total number of particles.
     * @return Raw pointer to acceleration vector buffer.
     */
    Tensor* get_accelerations(std::vector<Tensor> &acceleration_vector, const Particle* begin, unsigned int size);

  protected:
    std::vector<const Particle*> particle_ordering;                 /**< Array of particle pointers sorted by Morton code. */
    std::vector<std::vector<std::unique_ptr<NodeI>>> nodes_vector; /**< 2D grid storing unique pointers to tree nodes by depth level. */

    /**
     * @brief Initializes tree hierarchy iteratively using Morton encoding.
     * @param begin Start particle pointer.
     * @param end End particle pointer.
     */
    void init_tree(const Particle* begin, const Particle* end);

  public:
    /**
     * @brief Prints root node multipole coefficients to standard output.
     */
    void print_root_multipoles();

    /**
     * @brief Returns particle spatial partition depth and node ID tuples.
     * @return Vector of `(const Particle*, depth, node_id)`.
     */
    std::vector<std::tuple<const Particle*, int, int>> get_partition();

    /**
     * @brief Returns node interaction list tuples for visualization.
     * @return Vector of `(source_node, target_node, depth, interaction_type, particle)`.
     */
    std::vector<std::tuple<NodeI*, NodeI*, unsigned int, char, const Particle*>> get_nodes_interactions();
};