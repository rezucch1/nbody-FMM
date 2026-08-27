#pragma once

#include <memory>
#include "system_update_method.hpp"

/**
 * @file naive_update.hpp
 * @brief Direct \f$ \mathcal{O}(N^2) \f$ particle-particle gravitational acceleration solver.
 */

/**
 * @class NaiveUpdate
 * @brief Direct all-pairs \f$ \mathcal{O}(N^2) \f$ gravitational acceleration solver.
 *
 * Calculates exact gravitational acceleration via Newton's law of universal gravitation:
 * \f[ \mathbf{a}_i = G \sum_{j \neq i} m_j \frac{\mathbf{x}_j - \mathbf{x}_i}{\|\mathbf{x}_j - \mathbf{x}_i\|^3} \f]
 */
class NaiveUpdate : public SystemUpdateMethod {
  public:
    /**
     * @brief Computes exact \f$ \mathcal{O}(N^2) \f$ accelerations for all particles.
     * @param particles System particle vector.
     * @return Pointer to buffer of computed acceleration Tensors.
     */
    virtual const Tensor* update(const std::vector<Particle> &particles) override;

  protected:
    std::vector<Tensor> next_acceleration; /**< Buffer storing computed acceleration Tensors. */
};