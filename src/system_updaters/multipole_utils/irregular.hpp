#pragma once

#include "power_set.hpp"

/**
 * @file irregular.hpp
 * @brief Irregular solid harmonics \f$ I_l^m(\mathbf{r}) \f$ for 3D multipole conversions.
 */

/**
 * @class Irregular
 * @brief Class representing Irregular Solid Harmonics \f$ I_l^m(\mathbf{r}) \f$ used in 3D M2L operations.
 *
 * Evaluates singular/irregular spherical harmonics at displacement vector \f$ \mathbf{r} \f$:
 * \f[ I_l^m(\mathbf{r}) = \sqrt{\frac{4\pi}{2l+1}} \frac{Y_l^m(\theta, \phi)}{r^{l+1}} \f]
 */
class Irregular : public PowerSetI {
  public:
    /**
     * @brief Constructs Irregular solid harmonics set up to degree \f$ L \f$ for displacement vector \f$ \mathbf{r} \f$.
     * @param L Truncation degree limit.
     * @param r Displacement vector \f$ \mathbf{r} \f$.
     */
    Irregular(unsigned int L, const Tensor &r);

    /**
     * @brief Accesses complex Irregular Solid Harmonic coefficient \f$ I_l^m(\mathbf{r}) \f$.
     * @param l Harmonic degree \f$ l \f$.
     * @param m Harmonic order \f$ m \f$.
     * @return Complex harmonic value.
     */
    virtual std::complex<double> operator()(unsigned l, int m = 0) const override;
};
