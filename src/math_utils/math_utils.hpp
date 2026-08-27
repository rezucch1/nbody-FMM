#pragma once
#include <tuple>

/**
 * @file math_utils.hpp
 * @brief Utility functions for complex arithmetic operations.
 */

/**
 * @brief Performs inline complex multiplication of two numbers represented by real and imaginary components.
 * \f[ (a + ib)(c + id) = (ac - bd) + i(ad + bc) \f]
 * @param Re1 Real part of first complex number \f$ a \f$.
 * @param Im1 Imaginary part of first complex number \f$ b \f$.
 * @param Re2 Real part of second complex number \f$ c \f$.
 * @param Im2 Imaginary part of second complex number \f$ d \f$.
 * @return Tuple `{Real, Imaginary}` of resulting product.
 */
inline std::tuple<double, double> mult_complex(double Re1, double Im1, double Re2, double Im2) {
    return { Re1 * Re2 - Im1 * Im2, Re1 * Im2 + Re2 * Im1 };
}
