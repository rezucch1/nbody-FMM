#pragma once

#include <vector>

class SystemUpdateMethod;
#include "math_utils/tensor.hpp"
#include "particle.hpp"

/**
 * @file system_update_method.hpp
 * @brief Abstract base class for particle acceleration update solvers in gravitational simulations.
 */

/**
 * @brief Universal gravitational constant \f$ G = 6.6743 \times 10^{-11} \, \text{m}^3 \text{kg}^{-1} \text{s}^{-2} \f$.
 */
constexpr double GRAV_CONST = 6.6743e-11;

/**
 * @class SystemUpdateMethod
 * @brief Abstract interface for computing accelerations acting on particles.
 */
class SystemUpdateMethod {
  public:
    /**
     * @brief Computes acceleration vectors for a list of particles.
     * @param particles Const reference to vector of system particles.
     * @return Pointer to contiguous buffer/vector of computed acceleration Tensors.
     */
    virtual const Tensor* update(const std::vector<Particle> &particles) = 0;
};