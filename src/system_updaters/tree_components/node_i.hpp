#pragma once
#include <memory>
#include "power_set.hpp"
#include "particle.hpp"
class Node;
// #include "node.hpp"

//we need an Interface, with Multipole as a subclass (same with Leaf)
class NodeI{
  protected:
    NodeI(const Node *parent, Particle** particles_begin, Particle** particles_end)
    : parent(parent)
    , particles_begin(particles_begin)
    , particles_end(particles_end)
    , dim((*particles_begin)->get_position().dim)
    {}
    const Node *parent; //pointer to parent multipole (both multipole and leaf have a parent)

    Particle** const particles_begin;
    Particle** const particles_end;
    
    const unsigned int dim;

    Tensor mass_center;
    std::unique_ptr<MultipoleSetI> multipole_set;

  public:
    static const std::unique_ptr<MultipoleSetI> &getMultipoleSet(const NodeI &_this); 
    virtual void get_partition(std::vector<std::tuple<Particle*, int, int>> &partitions, int level, int partition_id) const;
  
  protected:
    static const Tensor &getMassCenter(const NodeI &_this); 
    
    virtual void calculateMC() = 0; //mass center calculator

};

