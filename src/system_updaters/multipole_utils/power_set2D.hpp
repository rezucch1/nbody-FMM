#pragma once
#include "power_set.hpp"

/**
 * @file power_set2D.hpp
 * @brief 2D PowerSet specialization for complex powers \f$ (z - z_0)^l \f$.
 */

/**
 * @class PowerSet<2>
 * @brief 2D specialization of PowerSet representing 2D complex coordinate powers.
 *
 * Computes powers of complex position offset \f$ z = x + i y \f$:
 * \f[ R_k(z) = z^k = (x + i y)^k \f]
 */
template <>
class PowerSet<2> : public PowerSetI {
  public:
    /**
     * @brief Constructs 2D PowerSet up to degree \f$ L \f$ for position displacement vector \f$ \mathbf{r} \f$.
     * @param L Truncation degree limit.
     * @param r Displacement vector \f$ \mathbf{r} = (x, y) \f$.
     */
    PowerSet(unsigned int L, const Tensor &r);

    /**
     * @brief Accesses complex power expansion coefficient.
     * @param l Power degree index \f$ l \f$.
     * @param m Unused in 2D (defaults to 0).
     * @return Complex value \f$ z^l \f$.
     */
    virtual std::complex<double> operator()(unsigned l, int m = 0) const override;
};

extern template class PowerSet<2>;
