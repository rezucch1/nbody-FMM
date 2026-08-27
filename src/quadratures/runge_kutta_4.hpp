#pragma once
#include "quadrature_method.hpp"

/**
 * @file runge_kutta_4.hpp
 * @brief Classical 4th-order Runge-Kutta (RK4) numerical integrator.
 */

/**
 * @class RungeKutta4
 * @brief Classical RK4 integrator for high-accuracy time integration of particle trajectories.
 *
 * Computes intermediate slopes \f$ k_1, k_2, k_3, k_4 \f$ to achieve 4th-order local truncation error.
 */
class RungeKutta4 : public QuadratureMethod {
  public:
    /**
     * @brief Default constructor.
     */
    RungeKutta4() {}

    /**
     * @brief Performs one 4th-order Runge-Kutta integration step.
     * @param method System acceleration update method pointer.
     * @param particles Particle vector to update.
     * @param delta_t Time step \f$ \Delta t \f$.
     */
    void integrate(SystemUpdateMethod *method, std::vector<Particle> &particles, const double &delta_t) override;
};