#pragma once

#include <memory>
#include <list>

#include "tensor.hpp"

class Multipole;

class Multipole{
  public:

  protected:
    Multipole *parent;
    std::unique_ptr<Multipole[]> children;
    
    Tensor mass_center;
    std::list<double[]> dipole_mometums;

};