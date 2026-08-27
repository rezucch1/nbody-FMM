#pragma once
#include "power_set.hpp"

/**
 * @file power_set3D.hpp
 * @brief 3D PowerSet specialization for Regular Solid Harmonics \f$ R_l^m(\mathbf{r}) \f$.
 */

/**
 * @class PowerSet<3>
 * @brief 3D specialization of PowerSet representing Regular Solid Harmonics \f$ R_l^m(\mathbf{r}) \f$.
 *
 * Computes regular solid harmonics evaluated at displacement vector \f$ \mathbf{r} = (x, y, z) \f$:
 * \f[ R_l^m(\mathbf{r}) = \sqrt{\frac{4\pi}{2l+1}} r^l Y_l^m(\theta, \phi) \f]
 */
template <>
class PowerSet<3> : public PowerSetI {
  public:
    /**
     * @brief Constructs 3D PowerSet up to degree \f$ L \f$ for position displacement vector \f$ \mathbf{r} \f$.
     * @param L Truncation degree limit.
     * @param r Displacement vector \f$ \mathbf{r} = (x, y, z) \f$.
     */
    PowerSet(unsigned int L, const Tensor &r);

    /**
     * @brief Accesses Regular Solid Harmonic coefficient \f$ R_l^m(\mathbf{r}) \f$.
     * @param l Harmonic degree \f$ l \f$.
     * @param m Harmonic order \f$ m \f$ (\f$ -l \le m \le l \f$).
     * @return Complex harmonic value.
     */
    virtual std::complex<double> operator()(unsigned l, int m = 0) const override;
};
