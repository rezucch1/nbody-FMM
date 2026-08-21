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
    Node(std::vector<std::vector<std::unique_ptr<NodeI>>> &allocator, unsigned int depth, unsigned int id_child, Particle** particles_begin, Particle** particles_end, const Tensor &a, const Tensor &b); //constructor
    // const std::vector<std::unique_ptr<NodeI>> &get_children() const;
    const std::vector<NodeI *> &get_neighbours() const;

    void get_partition(std::vector<std::tuple<Particle*, int, int>> &partitions) const override;

  protected:
    NodeI* get_child(unsigned int idx);
    std::vector<std::unique_ptr<NodeI>>::iterator get_children_begin();
    std::vector<std::unique_ptr<NodeI>>::iterator get_children_end();

    std::vector<std::unique_ptr<NodeI>>::const_iterator get_children_begin() const;
    std::vector<std::unique_ptr<NodeI>>::const_iterator get_children_end() const;
    virtual void calculateMC() override;
    //we're in the middle
};