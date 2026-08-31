#pragma once

#include "ghost_node_i.hpp"

/**
 * @file ghost_node.hpp
 * @brief Header for GhostNode representing remote internal tree nodes owned by other MPI ranks.
 */

/**
 * @class GhostNode
 * @brief Ghost node representing a remote internal tree cell.
 */
class GhostNode : public GhostNodeI {
  public:
    /**
     * @brief Constructs a GhostNode.
     * @param allocator Reference to global tree node grid allocator.
     * @param depth Tree depth level.
     * @param id_child Morton child ID.
     * @param dim Spatial dimension.
     */
    GhostNode(std::vector<std::unordered_map<uint64_t, std::unique_ptr<NodeI>>> &allocator, unsigned int depth, uint64_t id_child, unsigned int dim);

    virtual ~GhostNode() = default;

    /**
     * @brief No-op: Ghost nodes do not calculate M2M multipoles locally.
     * @param L Degree limit.
     */
    virtual void compute_multipoles(unsigned int L) override {};

    /**
     * @brief No-op: Ghost nodes do not collect M2L local expansions.
     */
    virtual void collect_multipoles_to_locals() override {};

    /**
     * @brief No-op: Ghost nodes do not propagate L2L local expansions.
     * @param parent_local Parent local expansion.
     */
    virtual void propagate_locals(const LocalSetI *parent_local = nullptr) override {};

  protected:
    /**
     * @brief Center of mass calculation for GhostNode.
     */
    virtual void calculateMC() override {};
};
