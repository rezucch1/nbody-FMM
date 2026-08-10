#pragma once
#include "math_utils/tensor.hpp"
class Particle;
#include "quadratures/quadrature_method.hpp"

class Particle{
  public:
    Particle(){};
    Particle(double mass, Tensor &&pos, Tensor &&vel);

    const double &get_mass() const;
    const Tensor &get_position() const;
    const Tensor &get_velocity() const;

    friend class QuadratureMethod;

  protected:
    double mass;
    Tensor position;
    Tensor velocity;
    
};