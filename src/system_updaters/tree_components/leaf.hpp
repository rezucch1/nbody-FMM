#pragma once
#include "node_i.hpp"
#include "particle.hpp"

class Leaf: public NodeI{
  protected:
    virtual void calculateMC() override;
  public:
    Leaf(const Node *parent, Particle** particles_begin,
         Particle** particles_end); //constructor
  




};