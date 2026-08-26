#pragma once
#include "math_utils/tensor.hpp"
class Particle;
#include "quadratures/quadrature_method.hpp"

class Particle{
  public:
    Particle(){};
    Particle(double weight, Tensor &&pos, Tensor &&vel);

    const Tensor get_acceleration(const Tensor &potential_gradiend) const;
    double get_weight() const;
    const Tensor &get_position() const;
    const Tensor &get_velocity() const;

    friend class QuadratureMethod;

  protected:
    double weight;
    Tensor position;
    Tensor velocity;
    
};