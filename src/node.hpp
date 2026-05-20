#pragma once //similar to "if not defined, define"

#include <memory>
#include <list>
#include "particle.hpp"
#include "tensor.hpp"

class Node; //because I need to use the class INSIDE the class
#include "multi_inter.hpp"


class Node: public Multi_inter{ //Multipole is a subclass of Multi_inter
  public:
    Node(const Node *parent, std::vector<Particle>::iterator &particles_begin, std::vector<Particle>::iterator &particles_end, const Tensor &a, const Tensor &b); //constructor

  protected:
    // const Multipole *parent; //pointer to parent multipole
    std::vector<std::unique_ptr<Multi_inter>> children; //array of children
    //we're in the middle
};