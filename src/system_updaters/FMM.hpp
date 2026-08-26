#pragma once

#include <memory>
#include "system_update_method.hpp"

class FMM : public SystemUpdateMethod{
  public:
    virtual const Tensor* update(const std::vector<Particle> &particles) override;

  protected:
    std::vector<Tensor> next_acceleration;
};