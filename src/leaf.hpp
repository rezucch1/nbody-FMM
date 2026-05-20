#pragma once
#include "multi_inter.hpp"

class Leaf: public Multi_inter{
  protected:
    const std::vector<Particle>::iterator particles_begin;
    const std::vector<Particle>::iterator particles_end;
    virtual void calculateMC() override;
  public:
    Leaf(const Node *parent, const std::vector<Particle>::iterator &particles_begin,
         const std::vector<Particle>::iterator &particles_end); //constructor
  




};