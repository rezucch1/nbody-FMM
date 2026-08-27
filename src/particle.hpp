#pragma once
#include "math_utils/tensor.hpp"
class Particle;
#include "quadratures/quadrature_method.hpp"

/**
 * @file particle.hpp
 * @brief Representation of physical particles in gravitational/field simulations.
 */

/**
 * @class Particle
 * @brief Class representing a point particle with mass/weight, position, and velocity vectors.
 */
class Particle {
  public:
    /**
     * @brief Default constructor for Particle.
     */
    Particle() {};

    /**
     * @brief Constructs a Particle with given weight/mass, position, and velocity.
     * @param weight Mass or weight parameter \$ m_i \$.
     * @param pos Position vector \$ \mathbf{x}_i \$.
     * @param vel Velocity vector \$ \mathbf{v}_i \$.
     */
    Particle(double weight, Tensor &&pos, Tensor &&vel);

    /**
     * @brief Calculates gravitational acceleration from potential gradient.
     * \$[ \mathbf{a}_i = -\nabla \Phi(\mathbf{x}_i) \$]
     * @param potential_gradiend Evaluated potential gradient vector \f$ \nabla \Phi \f$.
     * @return Acceleration vector \$ \mathbf{a}_i \$.
     */
    const Tensor get_acceleration(const Tensor &potential_gradiend) const;

    /**
     * @brief Gets particle mass/weight.
     * @return Mass or weight \f$ m_i \f$.
     */
    double get_weight() const;

    /**
     * @brief Gets reference to particle position vector.
     * @return Position vector \f$ \mathbf{x}_i \f$.
     */
    const Tensor &get_position() const;

    /**
     * @brief Gets reference to particle velocity vector.
     * @return Velocity vector \f$ \mathbf{v}_i \f$.
     */
    const Tensor &get_velocity() const;

    friend class QuadratureMethod;

  protected:
    double weight;   /**Particle mass or weight parameter. */
    Tensor position; /**Spatial position vector. */
    Tensor velocity; /**Spatial velocity vector. */
};