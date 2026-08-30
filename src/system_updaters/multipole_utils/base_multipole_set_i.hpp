#pragma once
#include <vector>
#include <complex>

/**
 * @file base_multipole_set_i.hpp
 * @brief Base class BaseMultipoleSetI for holding multipole and local expansion coefficient storage.
 */

/**
 * @class BaseMultipoleSetI
 * @brief Abstract base class holding expansion element coefficients up to degree \f$ L \f$.
 */
class BaseMultipoleSetI {
  protected:
    std::vector<double> elements; /**< Internal coefficient vector storing real/imaginary expansion components. */
    unsigned int L;               /**< Truncation degree limit \f$ L \f$ of expansion. */

    /**
     * @brief Constructs BaseMultipoleSetI with degree limit \f$ L \f$.
     * @param L Truncation degree.
     */
    BaseMultipoleSetI(unsigned int L);

  public:
    virtual ~BaseMultipoleSetI() = default;

    /**
     * @brief Accesses expansion element \f$ M_l^m \f$ or \f$ L_l^m \f$ at degree \f$ l \f$ and order \f$ m \f$.
     * @param l Degree index \f$ l \f$.
     * @param m Order index \f$ m \f$.
     * @return Complex expansion coefficient value.
     */
    virtual std::complex<double> operator()(unsigned l, int m = 0) const = 0;

    friend class MultipoleSetI;
};