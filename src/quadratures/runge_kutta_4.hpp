#pragma once
#include "quadrature_method.hpp"

class RungeKutta4 : public QuadratureMethod{
  public:
    RungeKutta4(){}
    void integrate(SystemUpdateMethod *method, std::vector<Particle> &particles, const double &delta_t) override;
};