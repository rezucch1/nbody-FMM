#pragma once

#include "local_set.hpp"

/**
 * @file local_set2D.hpp
 * @brief 2D LocalSet specialization representing 2D local Taylor expansions.
 */

/**
 * @class LocalSet<2>
 * @brief 2D specialization of LocalSet.
 *
 * Potential expansion:
 * \f[ \Phi(z) = \text{Re} \left( \sum_{l=0}^L b_l (z - z_0)^l \right) \f]
 * Gradient evaluation:
 * \f[ \nabla \Phi(z) = \left( \text{Re}\left( \sum_{l=1}^L l b_l z^{l-1} \right), -\text{Im}\left( \sum_{l=1}^L l b_l z^{l-1} \right) \right) \f]
 */
template <>
class LocalSet<2> : public LocalSetI {
  public:
    /**
     * @brief Constructs 2D LocalSet up to degree \f$ L \f$.
     * @param L Truncation degree limit.
     */
    LocalSet(unsigned int L);

    /**
     * @brief Gets complex local coefficient at degree \f$ l \f$.
     * @param l Degree \f$ l \f$.
     * @param m Unused in 2D (defaults to 0).
     * @return Complex local coefficient.
     */
    virtual std::complex<double> operator()(unsigned l, int m = 0) const override;

    /**
     * @brief Gets mutable reference to element at index \f$ i \f$.
     * @param i Coefficient index.
     * @return Mutable reference to complex coefficient.
     */
    std::complex<double> &get_element(unsigned int i);

    /**
     * @brief Gets const reference to element at index \f$ i \f$.
     * @param i Coefficient index.
     * @return Const reference to complex coefficient.
     */
    const std::complex<double> &get_element(unsigned int i) const;

    /**
     * @brief Polymorphic in-place addition of local expansion sets.
     * @param other Pointer to another LocalSetI.
     * @return Reference to this 2D LocalSet.
     */
    virtual LocalSet<2> &operator+=(const LocalSetI *other) override;

    /**
     * @brief In-place addition of 2D LocalSet.
     * @param other Target 2D LocalSet.
     * @return Reference to this 2D LocalSet.
     */
    virtual LocalSet<2> &operator+=(const LocalSet<2> other);

    /**
     * @brief Performs 2D L2L local shift translation to child cell center.
     * \f[ b_l^{\text{child}} = \sum_{k=l}^L \binom{k}{l} b_k^{\text{parent}} d^{k-l} \f]
     * @param d Displacement vector \f$ \mathbf{d} = \mathbf{x}_{\text{child}} - \mathbf{x}_{\text{parent}} \f$.
     * @return Pointer to translated 2D LocalSet.
     */
    virtual LocalSet<2> *distribute_parent_with_distance(const Tensor &d) const override;

    /**
     * @brief Evaluates 2D potential gradient \f$ \nabla \Phi \f$ at target offset \f$ \mathbf{d} \f$.
     * @param d Displacement vector from local expansion center to particle.
     * @return 2D gradient tensor.
     */
    virtual Tensor get_gradient(const Tensor &d) const override;

  protected:
    std::vector<std::complex<double>> elements; /**< Complex coefficient array storing \f$ b_l \f$. */
};

extern template class LocalSet<2>;
