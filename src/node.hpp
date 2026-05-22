#pragma once //similar to "if not defined, define"

#include <memory>
#include <list>
#include "particle.hpp"
#include "tensor.hpp"

class Node; //because I need to use the class INSIDE the class
#include "node_i.hpp"


class Node: public NodeI{ //Multipole is a subclass of NodeI
  public:
    Node(std::nullptr_t, __gnu_cxx::__normal_iterator<Particle *, std::vector<Particle, std::allocator<Particle>>>, __gnu_cxx::__normal_iterator<Particle *, std::vector<Particle, std::allocator<Particle>>>); //constructor
    Node(const Node *parent, std::vector<Particle>::iterator &particles_begin, std::vector<Particle>::iterator &particles_end, const Tensor &a, const Tensor &b); //constructor

  protected:
    void _constructor(std::vector<Particle>::iterator &particles_begin, std::vector<Particle>::iterator &particles_end, const Tensor &a, const Tensor &b);
    // const Multipole *parent; //pointer to parent multipole
    std::vector<std::unique_ptr<NodeI>> children; //array of children
    virtual void calculateMC() override;
    //we're in the middle
};