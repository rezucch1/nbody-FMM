#pragma once
#include "node_i.hpp"

class Leaf: public NodeI{
  protected:
    const std::vector<Particle>::iterator particles_begin;
    const std::vector<Particle>::iterator particles_end;
    virtual void calculateMC() override;
  public:
    Leaf(const Node *parent, const std::vector<Particle>::iterator &particles_begin,
         const std::vector<Particle>::iterator &particles_end); //constructor
  




};