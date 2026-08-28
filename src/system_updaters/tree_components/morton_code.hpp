#pragma once
#include <cstddef>

/**
 * @file morton_code.hpp
 * @brief Morton Z-order curve encoding and decoding helper utilities for multi-dimensional spatial indexing.
 */

/**
 * @class Morton
 * @brief Class for inter-leaving spatial coordinate bits to calculate Morton Z-order curve indices.
 */
class Morton {
  public:
    /**
     * @brief Encodes multi-dimensional integer Cartesian coordinates into a single 1D Morton code (Z-index).
     * @param cartesian_coordinates Array of integer grid coordinates \f$ [x_0, x_1, \dots, x_{d-1}] \f$.
     * @param dim Spatial dimension \f$ d \f$.
     * @return Encoded Morton Z-order integer code.
     */
    static size_t encode(const unsigned int* cartesian_coordinates, unsigned int dim);

    /**
     * @brief Decodes a 1D Morton code back into multi-dimensional integer Cartesian coordinates.
     * @param cartesian_coordinates Output array for decoded grid coordinates.
     * @param dim Spatial dimension \f$ d \f$.
     * @param code Morton Z-order integer code to decode.
     */
    static void decode(unsigned int* cartesian_coordinates, unsigned int dim, size_t code);
};