#pragma once
#include "local_set.hpp"

/**
 * @file local_set3D.hpp
 * @brief 3D LocalSet specialization representing 3D Spherical/Solid Harmonic local expansions.
 */

/**
 * @class LocalSet<3>
 * @brief 3D specialization of LocalSet.
 *
 * Manages 3D local expansion coefficients \f$ L_l^m \f$ and evaluates 3D potential gradients \f$ \nabla \Phi \f$.
 */
template <>
class LocalSet<3> : public LocalSetI {
  public:
    /**
     * @brief Constructs 3D LocalSet up to degree \f$ L \f$.
     * @param L Truncation degree limit.
     */
    LocalSet(unsigned int L);

    /**
     * @brief Creates a deep copy clone of the 3D local expansion set.
     * @return Unique pointer to cloned LocalSetI instance.
     */
    virtual std::unique_ptr<LocalSetI> clone() const override;

    /**
     * @brief Gets complex local coefficient \f$ L_l^m \f$.
     * @param l Degree \f$ l \f$.
     * @param m Order \f$ m \f$.
     * @return Complex coefficient.
     */
    virtual std::complex<double> operator()(unsigned l, int m = 0) const override;

    /**
     * @brief Polymorphic in-place addition of local expansion sets.
     * @param other Pointer to another LocalSetI.
     * @return Reference to this 3D LocalSet.
     */
    virtual LocalSet<3> &operator+=(const LocalSetI *other) override;

    /**
     * @brief In-place addition of 3D LocalSet.
     * @param other Target 3D LocalSet.
     * @return Reference to this 3D LocalSet.
     */
    virtual LocalSet<3> &operator+=(const LocalSet<3> other);

    /**
     * @brief Performs 3D L2L local shift translation to child cell center.
     * @param d Displacement vector \f$ \mathbf{d} = \mathbf{x}_{\text{child}} - \mathbf{x}_{\text{parent}} \f$.
     * @return Unique pointer to translated 3D LocalSet.
     */
    virtual std::unique_ptr<LocalSetI> distribute_parent_with_distance(const Tensor &d) const override;

    /**
     * @brief Sets element coefficient at degree \f$ n \f$ and order \f$ m \f$.
     * @param n Degree \f$ n \f$.
     * @param m Order \f$ m \f$.
     * @param val Complex value to set.
     */
    void set_elements(unsigned int n, int m, std::complex<double> val);

    /**
     * @brief Evaluates 3D potential gradient \f$ \nabla \Phi \f$ at target offset \f$ \mathbf{d} \f$.
     * @param d Displacement vector from local expansion center to evaluation point.
     * @return 3D gradient tensor.
     */
    virtual Tensor get_gradient(const Tensor &d) const override;
};

extern template class LocalSet<3>;