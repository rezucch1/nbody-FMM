#pragma once

#include <memory>
#include <list>
#include "particle.hpp"
#include "math_utils/tensor.hpp"

class Node;
#include "node_i.hpp"

/**
 * @file node.hpp
 * @brief Internal non-leaf Node class in the FMM spatial tree representation.
 */

/**
 * @class Node
 * @brief Represents an internal tree node managing child nodes and aggregate multipole expansions.
 */
class Node : public NodeI {
  public:
    /**
     * @brief Constructs an internal Node and initializes child nodes.
     * @param allocator Reference to 2D nodes grid.
     * @param depth Tree level depth.
     * @param id_child Morton child ID.
     * @param particles_begin Start particle pointer range.
     * @param particles_end End particle pointer range.
     */
    Node(std::vector<std::vector<std::unique_ptr<NodeI>>> &allocator, unsigned int depth, unsigned int id_child, const Particle** particles_begin, const Particle** particles_end);

    /**
     * @brief Computes multipole expansions recursively upward from children (M2M phase).
     * @param L Maximum degree of expansion.
     */
    virtual void compute_multipoles(unsigned int L) override;

    /**
     * @brief Collects multipoles from interaction list into local expansion set (M2L phase).
     */
    virtual void collect_multipoles_to_locals() override;

    /**
     * @brief Propagates parent local expansion down to children (L2L phase).
     * @param parent_local Translated local expansion from parent node.
     */
    virtual void propagate_locals(LocalSetI *parent_local = nullptr) override;

    /**
     * @brief Gets list of neighbor nodes.
     * @return Reference to neighbor node pointers vector.
     */
    const std::vector<NodeI *> &get_neighbours() const;

    /**
     * @brief Recursively collects particle spatial partition metadata.
     * @param partitions Reference to partition output tuple list.
     */
    void get_partition(std::vector<std::tuple<const Particle*, int, int>> &partitions) const override;

  protected:
    /**
     * @brief Gets child node pointer by index.
     * @param idx Child index (0 to \f$ 2^d - 1 \f$).
     * @return Pointer to child node.
     */
    NodeI* get_child(unsigned int idx);

    /**
     * @brief Gets iterator to start of child nodes.
     * @return Iterator.
     */
    std::vector<std::unique_ptr<NodeI>>::iterator get_children_begin();

    /**
     * @brief Gets iterator to end of child nodes.
     * @return Iterator.
     */
    std::vector<std::unique_ptr<NodeI>>::iterator get_children_end();

    /**
     * @brief Gets const iterator to start of child nodes.
     * @return Const iterator.
     */
    std::vector<std::unique_ptr<NodeI>>::const_iterator get_children_begin() const;

    /**
     * @brief Gets const iterator to end of child nodes.
     * @return Const iterator.
     */
    std::vector<std::unique_ptr<NodeI>>::const_iterator get_children_end() const;

    /**
     * @brief Calculates center of mass from child node centers of mass.
     */
    virtual void calculateMC() override;
};