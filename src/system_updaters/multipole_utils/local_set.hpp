#pragma once
#include <memory>
#include "base_multipole_set_i.hpp"
#include "math_utils/tensor.hpp"

/**
 * @file local_set.hpp
 * @brief Base LocalSetI class and templates for 2D/3D local expansions \f$ L_l^m \f$.
 */

/**
 * @class LocalSetI
 * @brief Abstract interface representing local expansions (Taylor expansions around cell center) and evaluation of potential gradients.
 */
class LocalSetI : public BaseMultipoleSetI {
  public:
    /**
     * @brief Creates a deep copy clone of the local expansion set.
     * @return Unique pointer to cloned LocalSetI instance.
     */
    virtual std::unique_ptr<LocalSetI> clone() const = 0;

    /**
     * @brief Polymorphic in-place addition of local expansion sets.
     * @param other Pointer to another LocalSetI.
     * @return Reference to updated LocalSetI.
     */
    virtual LocalSetI &operator+=(const LocalSetI *other) = 0;

    /**
     * @brief Translates parent local expansion center to child cell center (L2L step).
     * \f[ L_{\text{child}} = \text{L2L}(L_{\text{parent}}, \mathbf{d}) \f]
     * @param d Displacement vector \f$ \mathbf{d} = \mathbf{x}_{\text{child}} - \mathbf{x}_{\text{parent}} \f$.
     * @return Unique pointer to dynamically allocated translated LocalSetI.
     */
    virtual std::unique_ptr<LocalSetI> distribute_parent_with_distance(const Tensor &d) const = 0;

    /**
     * @brief Evaluates potential gradient \f$ \nabla \Phi(\mathbf{d}) \f$ at target offset \f$ \mathbf{d} \f$.
     * \f[ \mathbf{F} = -\nabla \Phi(\mathbf{d}) \f]
     * @param d Displacement vector from local expansion center to target evaluation point.
     * @return Gradient tensor.
     */
    virtual Tensor get_gradient(const Tensor &d) const = 0;

  protected:
    /**
     * @brief Constructs LocalSetI with degree limit \f$ L \f$.
     * @param L Truncation degree limit.
     */
    LocalSetI(unsigned int L) : BaseMultipoleSetI(L) {};
};

/**
 * @class LocalSet
 * @brief Primary template for dimension-specific LocalSet classes.
 * @tparam dim Spatial dimension (2 or 3).
 */
template <unsigned int dim>
class LocalSet : public LocalSetI {
};

#include "local_set2D.hpp"
#include "local_set3D.hpp"