#pragma once
#include "quadrature_method.hpp"

class SyplecticEuler : public QuadratureMethod{

  public:
    SyplecticEuler(){};
    void integrate(SystemUpdateMethod *method, std::vector<Particle> &particles, const double &delta_t) override;
  
};