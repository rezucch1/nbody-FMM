#pragma once
#include "node_i.hpp"
#include "particle.hpp"

class Leaf: public NodeI{
  protected:
    virtual void calculateMC() override;
  public:
    Leaf(std::vector<std::vector<std::unique_ptr<NodeI>>> &allocator, unsigned int depth, unsigned int child_id, const Particle** particles_begin,
         const Particle** particles_end); //constructor
    
    virtual void compute_multipoles(unsigned int L) override;

    virtual void propagate_locals(LocalSetI *parent_local = nullptr) override{
      NodeI::propagate_locals(parent_local);
    };

    void get_acceleration_vector(std::vector<Tensor> &acceleration);
  
};