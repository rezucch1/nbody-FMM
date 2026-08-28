#pragma once
#include <vector>
#include "math_utils/tensor.hpp"
#include <complex>
#include "base_multipole_set_i.hpp"
#include <iostream>
#include <memory>

class MultipoleSetI;
#include "power_set.hpp"
#include "local_set.hpp"

/**
 * @file multipole_set.hpp
 * @brief Base MultipoleSetI class and templates for 2D/3D multipole expansions \f$ M_l^m \f$.
 */

/**
 * @class MultipoleSetI
 * @brief Abstract class managing multipole moments and expansion transformations (M2M, M2L).
 */
class MultipoleSetI : public BaseMultipoleSetI {
  protected:
    /**
     * @brief Constructs MultipoleSetI with degree limit \f$ L \f$.
     * @param L Maximum degree.
     */
    MultipoleSetI(unsigned int L);

    friend std::ostream &operator<<(std::ostream &os, const MultipoleSetI &multipole);

  public:
    /**
     * @brief Accumulates power set terms (P2M step).
     * @param z PowerSet containing particle displacement powers.
     * @return Reference to updated MultipoleSetI.
     */
    virtual MultipoleSetI &operator+=(const PowerSetI &z);

    /**
     * @brief In-place addition of another multipole expansion.
     * @param z MultipoleSet to add.
     * @return Reference to updated MultipoleSetI.
     */
    MultipoleSetI &operator+=(const MultipoleSetI &z);

    /**
     * @brief Assigns null state (zeroes out expansion elements).
     * @return Reference to updated MultipoleSetI.
     */
    MultipoleSetI &operator=(std::nullptr_t);

    /**
     * @brief Translates child multipole expansion to parent center by distance \f$ \mathbf{d} \f$ (M2M step).
     * \f[ M_{\text{parent}} = \text{M2M}(M_{\text{child}}, \mathbf{d}) \f]
     * @param d Distance displacement vector from child center to parent center.
     * @return Unique pointer to translated multipole set.
     */
    virtual std::unique_ptr<MultipoleSetI> weigh_children_with_distance(const Tensor &d) const = 0;

    /**
     * @brief Converts multipole expansion into local expansion shifted by vector \f$ \mathbf{d} \f$ (M2L step).
     * \f[ L_{\text{local}} = \text{M2L}(M_{\text{source}}, \mathbf{d}) \f]
     * @param d Distance displacement vector between cell centers.
     * @return Unique pointer to converted local set.
     */
    virtual std::unique_ptr<LocalSetI> to_local(const Tensor &d) const = 0;

    /**
     * @brief Gets raw reference to internal elements array.
     * @param _this Target multipole set.
     * @return Elements vector reference.
     */
    static const std::vector<double> &getElements(const MultipoleSetI &_this);
};

/**
 * @class MultipoleSet
 * @brief Primary template for dimension-specific MultipoleSet classes.
 * @tparam dim Spatial dimension (2 or 3).
 */
template <unsigned int dim>
class MultipoleSet : public MultipoleSetI {
};

#include "multipole_set2D.hpp"
#include "multipole_set3D.hpp"
