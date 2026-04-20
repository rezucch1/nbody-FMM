#pragma once
#include "tensor.hpp"
class Particle;
#include "quadrature_method.hpp"

class Particle{
  public:
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