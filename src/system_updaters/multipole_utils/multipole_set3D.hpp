#pragma once

#include "multipole_set.hpp"

/**
 * @file multipole_set3D.hpp
 * @brief 3D MultipoleSet specialization using Spherical/Solid Harmonics \f$ M_l^m \f$.
 */

/**
 * @class MultipoleSet<3>
 * @brief 3D specialization of MultipoleSet.
 *
 * Manages 3D spherical multipole moments \f$ M_l^m \f$:
 * \f[ M_l^m = \sum_i m_i R_l^{-m}(\mathbf{r}_i) \f]
 */
template <>
class MultipoleSet<3> : public MultipoleSetI {
  public:
    /**
     * @brief Constructs 3D MultipoleSet initialized to zero.
     * @param L Truncation degree limit \f$ L \f$.
     */
    MultipoleSet(unsigned int L);

    /**
     * @brief Gets complex multipole moment \f$ M_l^m \f$.
     * @param l Harmonic degree \f$ l \f$.
     * @param m Harmonic order \f$ m \f$.
     * @return Complex multipole moment coefficient.
     */
    virtual std::complex<double> operator()(unsigned l, int m) const override;

    /**
     * @brief Accumulates particle 3D regular solid harmonic power terms (P2M step).
     * @param z PowerSet containing particle displacement harmonics.
     * @return Reference to updated 3D MultipoleSet.
     */
    virtual MultipoleSetI &operator+=(const PowerSetI &z) override;

    /**
     * @brief Performs 3D M2M multipole translation to parent center.
     * @param d Displacement vector \f$ \mathbf{d} = \mathbf{x}_{\text{parent}} - \mathbf{x}_{\text{child}} \f$.
     * @return Unique pointer to translated 3D MultipoleSet.
     */
    virtual std::unique_ptr<MultipoleSetI> weigh_children_with_distance(const Tensor &d) const override;

    /**
     * @brief Performs 3D M2L conversion using irregular solid harmonics.
     * @param d Displacement vector \f$ \mathbf{d} = \mathbf{x}_{\text{target}} - \mathbf{x}_{\text{source}} \f$.
     * @return Unique pointer to converted 3D LocalSet.
     */
    virtual std::unique_ptr<LocalSetI> to_local(const Tensor &d) const override;

    using MultipoleSetI::operator=;
};
