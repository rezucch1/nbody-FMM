#pragma once
#include <vector>
#include "math_utils/tensor.hpp"
#include <complex>
#include "base_multipole_set_i.hpp"

class PowerSetI;
#include "multipole_set.hpp"

/**
 * @file power_set.hpp
 * @brief PowerSet classes representing coordinate power terms and regular solid harmonics \f$ R_l^m(r) \f$.
 */

/**
 * @class PowerSetI
 * @brief Interface for particle coordinate power expansions used in P2M operations.
 */
class PowerSetI : public BaseMultipoleSetI {
  protected:
    /**
     * @brief Constructs PowerSetI with expansion degree \f$ L \f$.
     * @param L Truncation degree.
     */
    PowerSetI(unsigned int L);

  public:
    /**
     * @brief In-place weight/mass scaling of power set elements.
     * \f[ R_l^m \leftarrow m_i R_l^m \f]
     * @param weight Particle weight or mass \f$ m_i \f$.
     * @return Reference to this PowerSetI instance.
     */
    PowerSetI &operator*=(double weight);
};

/**
 * @class PowerSet
 * @brief Primary template for dimension-specific PowerSet implementations.
 * @tparam dim Spatial dimension (2 or 3).
 */
template <unsigned int dim>
class PowerSet : public PowerSetI {
};

#include "power_set2D.hpp"
#include "power_set3D.hpp"
