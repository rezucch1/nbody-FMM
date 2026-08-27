/**
 * @file quadrature_method.cpp
 * @brief Implementation of QuadratureMethod helper accessors.
 */

#include "quadrature_method.hpp"

Tensor &QuadratureMethod::position(Particle &p) { return p.position; }
Tensor &QuadratureMethod::velocity(Particle &p) { return p.velocity; }