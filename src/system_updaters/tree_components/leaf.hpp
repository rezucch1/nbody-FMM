#pragma once
#include "node_i.hpp"
#include "particle.hpp"

class Leaf: public NodeI{
  protected:
    virtual void calculateMC() override;
  public:
    Leaf(std::vector<std::vector<std::unique_ptr<NodeI>>> &allocator, unsigned int depth, unsigned int child_id, Particle** particles_begin,
         Particle** particles_end); //constructor
  




};