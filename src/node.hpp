#pragma once

#include <memory>
#include <list>
#include "particle.hpp"
#include "tensor.hpp"

class Multipole;

class Multipole{
  public:
    Multipole(const Multipole *parent, std::vector<Particle>::iterator &particles_begin, std::vector<Particle>::iterator &particles_end, const Tensor &a, const Tensor &b); //constructor

  protected:
    const Multipole *parent; //pointer to parent multipole
    std::vector<Multipole> children; //array of children
    //we're in the middle
    Tensor mass_center;
    std::list<double[]> dipole_mometums;

};