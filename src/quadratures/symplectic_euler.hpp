#pragma once
#include "quadrature_method.hpp"

/**
 * @file symplectic_euler.hpp
 * @brief Symplectic Euler 1st-order numerical integrator for Hamiltonian particle systems.
 */

/**
 * @class SyplecticEuler
 * @brief Symplectic Euler integration scheme preserving phase space geometric structure.
 *
 * Updates velocity first using current acceleration, then advances position using updated velocity:
 * \f[ \mathbf{v}_{n+1} = \mathbf{v}_n + \mathbf{a}_n \Delta t \f]
 * \f[ \mathbf{x}_{n+1} = \mathbf{x}_n + \mathbf{v}_{n+1} \Delta t \f]
 */
class SyplecticEuler : public QuadratureMethod {
  public:
    /**
     * @brief Default constructor.
     */
    SyplecticEuler() {};

    /**
     * @brief Performs one Symplectic Euler time step update.
     * @param method System acceleration update method pointer.
     * @param particles Particle vector to update in-place.
     * @param delta_t Time step \f$ \Delta t \f$.
     */
    void integrate(SystemUpdateMethod *method, std::vector<Particle> &particles, const double &delta_t) override;
};