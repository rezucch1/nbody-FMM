#pragma once
#include "node_i.hpp"
#include "particle.hpp"

/**
 * @file leaf.hpp
 * @brief Terminal Leaf node class representing spatial cells containing particle elements in FMM tree.
 */

/**
 * @class Leaf
 * @brief Terminal node in spatial tree that directly evaluates particle accelerations.
 */
class Leaf : public NodeI {
  protected:
    /**
     * @brief Calculates center of mass of particles directly contained in this leaf cell.
     * \f[ \mathbf{x}_{cm} = \frac{\sum_i m_i \mathbf{x}_i}{\sum_i m_i} \f]
     */
    virtual void calculateMC() override;
    const std::vector<const Particle *> particles;

  public:
    /**
     * @brief Constructs a Leaf node managing a specific particle pointer slice.
     * @param allocator Grid allocator reference.
     * @param depth Leaf level depth.
     * @param child_id Morton child index.
     * @param particles_begin Start particle pointer range.
     * @param particles_end End particle pointer range.
     */
    Leaf(std::vector<std::unordered_map<uint64_t, std::unique_ptr<NodeI>>> &allocator, unsigned int depth, uint64_t child_id, std::vector<const Particle *>::const_iterator particles_begin, std::vector<const Particle *>::const_iterator particles_end);

    /**
     * @brief Computes multipole expansion from particles in this leaf cell (P2M phase).
     * @param L Expansion degree limit.
     */
    virtual void compute_multipoles(unsigned int L) override;

    /**
     * @brief Propagates parent local expansion down to leaf.
     * @param parent_local Translated local expansion from parent node.
     */
    virtual void propagate_locals(const LocalSetI *parent_local = nullptr) override {
      NodeI::propagate_locals(parent_local);
    };

    /**
     * @brief Gets contained particles vector reference.
     * @return Const reference to particles array.
     */
    const std::vector<const Particle *> &get_particles() const { return particles; }

    /**
     * @brief Evaluates total acceleration vectors for all particles in this leaf cell.
     *
     * Combines far-field local potential expansion gradient with near-field direct particle-particle interactions from neighbor nodes.
     * @param acceleration Output acceleration vector list.
     */
    void compute_acceleration();
};