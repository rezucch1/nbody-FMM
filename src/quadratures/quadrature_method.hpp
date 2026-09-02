#pragma once

class QuadratureMethod;
#include "particle.hpp"
#include "system_updaters/system_update_method.hpp"

/**
 * @file quadrature_method.hpp
 * @brief Abstract base class for ODE numerical integration methods (quadratures).
 */

/**
 * @class QuadratureMethod
 * @brief Abstract interface for numerical integrators that advance particle states in time.
 */
class QuadratureMethod {
  public:
    /**
     * @brief Integrates particle positions and velocities over a time step \f$ \Delta t \f$.
     * @param method Pointer to SystemUpdateMethod for acceleration computation.
     * @param particles Vector of particles to advance in time.
     * @param delta_t Time step duration \f$ \Delta t \f$.
     */
    virtual void integrate(SystemUpdateMethod *method, std::vector<Particle> &particles, const double &delta_t) = 0;

  protected:
    /**
     * @brief Accesses mutable position tensor of a particle.
     * @param p Reference to particle.
     * @return Mutable reference to particle position.
     */
    Tensor &position(Particle &p);

    /**
     * @brief Accesses mutable velocity tensor of a particle.
     * @param p Reference to particle.
     * @return Mutable reference to particle velocity.
     */
    Tensor &velocity(Particle &p);

    /**
     * @brief Accesses mutable acceleration tensor of a particle.
     * @param p Reference to particle.
     * @return Mutable reference to particle acceleration.
     */
    Tensor &acceleration(Particle &p);
};