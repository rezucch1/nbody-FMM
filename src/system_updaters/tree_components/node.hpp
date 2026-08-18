#pragma once //similar to "if not defined, define"

#include <memory>
#include <list>
#include "particle.hpp"
#include "tensor.hpp"

class Node; //because I need to use the class INSIDE the class
#include "node_i.hpp"


class Node: public NodeI{ //Multipole is a subclass of NodeI
  public:
    // Node(const Node *parent, std::vector<Particle>::iterator &&particles_begin, std::vector<Particle>::iterator &&particles_end); //constructor

    // const Multipole *parent; //pointer to parent multipole
    Node(const Node *parent, Particle** particles_begin, Particle** particles_end, const Tensor &a, const Tensor &b); //constructor

    void get_partition(std::vector<std::tuple<Particle*, int, int>> &partitions, int level, int partition_id) const override;

  protected:
    std::vector<std::unique_ptr<NodeI>> children; //array of children
    virtual void calculateMC() override;
    //we're in the middle
};