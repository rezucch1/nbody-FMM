#pragma once

#include <vector>

class SystemUpdateMethod;
#include "math_utils/tensor.hpp"
#include "particle.hpp"

constexpr double GRAV_CONST = 6.6743e-11;

class SystemUpdateMethod{
  public:
    virtual const Tensor* update(const std::vector<Particle> &particles) = 0;
};