#pragma once

#include "node_i.hpp"

/**
 * @file ghost_node_i.hpp
 * @brief Abstract interface GhostNodeI for remote ghost tree cells (GhostNode and GhostLeaf).
 */

/**
 * @class GhostNodeI
 * @brief Interface for remote ghost cells providing set_multipole_set without needing type checks.
 */
class GhostNodeI : public NodeI {
  protected:
    /**
     * @brief Constructs GhostNodeI base instance.
     * @param allocator Reference to global tree node grid allocator.
     * @param depth Tree depth level.
     * @param id_child Morton child ID.
     * @param dim Spatial dimension.
     */
    GhostNodeI(std::vector<std::unordered_map<uint64_t, std::unique_ptr<NodeI>>> &allocator, unsigned int depth, uint64_t id_child, unsigned int dim)
    : NodeI(allocator, depth, id_child, dim) {}

  public:
    virtual ~GhostNodeI() = default;

    /**
     * @brief Identifies cell as ghost node.
     * @return Always true.
     */
    virtual bool is_ghost() const override { return true; }

    /**
     * @brief Sets remote multipole expansion set received via MPI.
     * @param ms Unique pointer to received MultipoleSetI.
     */
    virtual void set_multipole_set(std::unique_ptr<MultipoleSetI> ms) {
      multipole_set = std::move(ms);
      std::cout << "Multipole move with sucess\n";
    }
};
