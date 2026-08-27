#pragma once

#include "multipole_set.hpp"

/**
 * @file multipole_set2D.hpp
 * @brief 2D MultipoleSet specialization for 2D logarithmic/complex multipole expansions.
 */

/**
 * @class MultipoleSet<2>
 * @brief 2D specialization of MultipoleSet.
 *
 * Implements 2D M2M shift using binomial expansion:
 * \f[ M_l = \sum_{k=0}^{l} \binom{l}{k} d^k M_{l-k}^{\text{child}} \f]
 * and 2D M2L conversion to local Taylor expansion:
 * \f[ L_l = (-1)^l \sum_{k=0}^{L-l} \binom{l+k}{k} \frac{M_k}{d^{l+k+1}} \f]
 */
template <>
class MultipoleSet<2> : public MultipoleSetI {
  public:
    /**
     * @brief Constructs 2D MultipoleSet initialized to zero.
     * @param L Truncation degree limit \f$ L \f$.
     */
    explicit MultipoleSet(unsigned int L);

    /**
     * @brief Gets complex multipole coefficient at degree \f$ l \f$.
     * @param l Degree \f$ l \f$.
     * @param m Unused in 2D (defaults to 0).
     * @return Complex multipole moment.
     */
    virtual std::complex<double> operator()(unsigned l, int m = 0) const override;

    /**
     * @brief Performs 2D M2M multipole shift to parent cell center.
     * @param d Displacement vector \f$ \mathbf{d} = \mathbf{x}_{\text{parent}} - \mathbf{x}_{\text{child}} \f$.
     * @return Unique pointer to translated 2D MultipoleSet.
     */
    virtual std::unique_ptr<MultipoleSetI> weigh_children_with_distance(const Tensor &d) const override;

    /**
     * @brief Performs 2D M2L conversion to local expansion.
     * @param d Displacement vector between cell centers \f$ \mathbf{d} = \mathbf{x}_{\text{target}} - \mathbf{x}_{\text{source}} \f$.
     * @return Unique pointer to converted 2D LocalSet.
     */
    virtual std::unique_ptr<LocalSetI> to_local(const Tensor &d) const override;

    using MultipoleSetI::operator=;
};

extern template class MultipoleSet<2>;
