#pragma once
#include <memory>
#include <unordered_map>
#include <cstdint>
#include "system_updaters/multipole_utils/power_set.hpp"
#include "particle.hpp"
class Node;

/**
 * @file node_i.hpp
 * @brief Abstract base class NodeI representing an internal or leaf cell node in the FMM spatial tree hierarchy.
 */

/**
 * @class NodeI
 * @brief Abstract interface for spatial tree cells managing particle ranges, multipole sets, and local expansion sets.
 */
class NodeI {
  protected:
    /**
     * @brief Constructs NodeI base instance with depth, node ID, and particle pointer range.
     * @param allocator Reference to 2D nodes storage grid.
     * @param depth Depth level of node in tree hierarchy.
     * @param id_child Morton curve child index.
     * @param particles_begin Start pointer of contained particle array.
     * @param particles_end End pointer of contained particle array.
     */
    NodeI(std::vector<std::unordered_map<uint64_t, std::unique_ptr<NodeI>>> &allocator, const unsigned int depth, uint64_t id_child, unsigned int dim)
    : allocator(allocator)
    , depth(depth)
    , id_child(id_child)
    , dim(dim)
    {}

    std::vector<std::unordered_map<uint64_t, std::unique_ptr<NodeI>>> &allocator; /**< Reference to global tree node grid allocator. */
    const unsigned int depth;                                   /**< Tree depth level of this node. */
    const uint64_t id_child;                                /**< Morton code ID of this node at current depth. */

    std::vector<NodeI*> neighbours_list;  /**< List of adjacent spatial neighbor nodes. */
    std::vector<NodeI*> interaction_list; /**< List of well-separated interaction nodes (well-separated parent-neighbor children). */

    const unsigned int dim; /**< Spatial dimension (2 or 3). */

    Tensor mass_center;                            /**< Weighted center of mass position \f$ \mathbf{x}_{cm} \f$. */
    std::unique_ptr<MultipoleSetI> multipole_set;  /**< Multipole expansion moments set \f$ M_l^m \f$. */
    std::unique_ptr<LocalSetI> local_set;          /**< Local expansion coefficients set \f$ L_l^m \f$. */

  public:
    virtual ~NodeI() = default;

    /**
     * @brief Checks if this node is a remote ghost cell.
     * @return True if ghost node, false if local node.
     */
    virtual bool is_ghost() const { return false; }

    /**
     * @brief Retrieves parent Node reference.
     * @return Parent Node reference.
     */
    Node &get_parent();

    /**
     * @brief Static getter for a node's MultipoleSet reference.
     * @param _this Target node.
     * @return Reference to unique_ptr of MultipoleSetI.
     */
    static const std::unique_ptr<MultipoleSetI> &getMultipoleSet(const NodeI &_this); 

    /**
     * @brief Computes neighbor and interaction lists using Morton indexing arithmetic.
     */
    void compute_interaction_list();

    /**
     * @brief Gets the Morton node ID.
     * @return Morton child index.
     */
    unsigned int get_id() const;

    /**
     * @brief Pure virtual method to calculate multipole moments up to degree \f$ L \f$.
     * @param L Maximum expansion degree.
     */
    virtual void compute_multipoles(unsigned int L) = 0;

    /**
     * @brief Converts multipole expansions from interaction list nodes into local expansions (M2L phase).
     * \f[ L_{\text{local}} += \sum_{s \in \text{InteractionList}} \text{M2L}(M_s) \f]
     */
    virtual void collect_multipoles_to_locals();

    /**
     * @brief Propagates parent local expansion down to this node (L2L phase).
     * @param parent_local Translated local expansion from parent node.
     */
    virtual void propagate_locals(const LocalSetI *parent_local = nullptr);

  protected:
    /**
     * @brief Gets center of mass tensor reference from a node.
     * @param _this Target node.
     * @return Center of mass position tensor.
     */
    static const Tensor &getMassCenter(const NodeI &_this); 

    /**
     * @brief Pure virtual calculation of center of mass / weight center.
     */
    virtual void calculateMC() = 0;

    friend class Tree;
};
