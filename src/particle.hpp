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
 * @brief Class representing a point particle with mass/weight, position, velocity, and acceleration vectors.
 */
class Particle {
  public:
    /**
     * @brief Default constructor for Particle.
     */
    Particle() {};

    /**
     * @brief Constructs a Particle with given weight/mass, position, and velocity.
     * @param weight Mass or weight parameter \f$ m_i \f$.
     * @param pos Position vector \f$ \mathbf{x}_i \f$.
     * @param vel Velocity vector \f$ \mathbf{v}_i \f$.
     */
    Particle(double weight, Tensor &&pos, Tensor &&vel);

    /**
     * @brief Calculates and updates gravitational acceleration from potential gradient.
     * \f[ \mathbf{a}_i = - G \nabla \Phi(\mathbf{x}_i) \f]
     * @param potential_gradiend Evaluated potential gradient vector \f$ \nabla \Phi \f$.
     * @return Const reference to updated acceleration vector \f$ \mathbf{a}_i \f$.
     */
    const Tensor &compute_new_accelaration(const Tensor &potential_gradiend) const;

    /**
     * @brief Gets const reference to particle acceleration vector.
     * @return Acceleration vector \f$ \mathbf{a}_i \f$.
     */
    const Tensor &get_acceleration() const;

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
    double weight;               /**< Particle mass or weight parameter. */
    Tensor position;             /**< Spatial position vector. */
    Tensor velocity;             /**< Spatial velocity vector. */
    mutable Tensor acceleration; /**< Particle acceleration vector. */
};