#pragma once

#include "ghost_node_i.hpp"

/**
 * @file ghost_leaf.hpp
 * @brief Header for GhostLeaf representing remote leaf cells containing ghost particles from adjacent MPI ranks.
 */

/**
 * @class GhostLeaf
 * @brief Ghost leaf cell holding remote ghost particles.
 */
class GhostLeaf : public GhostNodeI {
  public:
    /**
     * @brief Constructs a GhostLeaf for a range of ghost particle pointers.
     * @tparam InputIt Iterator over Particle pointers.
     * @param allocator Reference to global tree node grid allocator.
     * @param depth Tree depth level.
     * @param id_child Morton child ID.
     * @param begin Start iterator over ghost particles.
     * @param end End iterator over ghost particles.
     */
    template <typename InputIt>
    GhostLeaf(std::vector<std::unordered_map<uint64_t, std::unique_ptr<NodeI>>> &allocator, unsigned int depth, uint64_t id_child, InputIt begin, InputIt end)
    : GhostNodeI(allocator, depth, id_child, begin < end ? (*begin)->get_position().dim : 3)
    {
        particles.reserve(end - begin);
        for (auto p = begin; p < end; ++p) {
            particles.push_back(*p);
        }
        calculateMC();
    }

    virtual ~GhostLeaf() = default;

    /**
     * @brief No-op: Ghost leaves do not compute P2M/M2M multipoles.
     * @param L Maximum expansion degree.
     */
    virtual void compute_multipoles(unsigned int L) override {};

    /**
     * @brief No-op: Ghost leaves do not collect M2L local expansions.
     */
    virtual void collect_multipoles_to_locals() override {};

    /**
     * @brief No-op: Ghost leaves do not propagate L2L local expansions.
     * @param parent_local Parent local expansion.
     */
    virtual void propagate_locals(const LocalSetI *parent_local = nullptr) override {};

    /**
     * @brief No-op: Ghost particles do not compute acceleration acting on themselves.
     */
    void compute_acceleration() {};

    /**
     * @brief Gets contained ghost particles vector reference.
     * @return Const reference to ghost particles array.
     */
    const std::vector<const Particle *> &get_particles() const { return particles; }

  protected:
    /**
     * @brief Calculates center of mass of contained ghost particles.
     */
    virtual void calculateMC() override;

  private:
    std::vector<const Particle *> particles; /**< Pointers to contained ghost particles. */
};
