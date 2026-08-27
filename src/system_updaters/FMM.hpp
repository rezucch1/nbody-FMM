#pragma once

#include <memory>
#include "system_update_method.hpp"

/**
 * @file FMM.hpp
 * @brief Fast Multipole Method (FMM) \f$ \mathcal{O}(N) / \mathcal{O}(N \log N) \f$ system updater.
 */

/**
 * @class FMM
 * @brief SystemUpdateMethod solver using the Fast Multipole Method (FMM).
 *
 * Constructs a spatial decomposition tree, evaluates multipole expansions, converts multipoles to local expansions,
 * and calculates potential gradients to achieve fast gravitational field updates.
 */
class FMM : public SystemUpdateMethod {
  public:
    /**
     * @brief Computes gravitational accelerations using Fast Multipole Method.
     * @param particles System particle vector.
     * @return Pointer to buffer of computed acceleration Tensors.
     */
    virtual const Tensor* update(const std::vector<Particle> &particles) override;

  protected:
    std::vector<Tensor> next_acceleration; /**< Buffer storing computed acceleration Tensors. */
};