#pragma once

class QuadratureMethod;
#include "particle.hpp"

#include "system_update_method.hpp"

class QuadratureMethod{
  public:
    virtual void integrate(SystemUpdateMethod *method, std::vector<Particle> &particles, const double &delta_t) = 0;
  protected:
    Tensor &position(Particle &p);
    Tensor &velocity(Particle &p);
};