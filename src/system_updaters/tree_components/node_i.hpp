#pragma once
#include <memory>
#include "power_set.hpp"
#include "particle.hpp"
class Node;
// #include "node.hpp"

//we need an Interface, with Multipole as a subclass (same with Leaf)
class NodeI{
  protected:
    NodeI(std::vector<std::vector<std::unique_ptr<NodeI>>> &allocator, const unsigned int depth, unsigned int id_child, Particle** particles_begin, Particle** particles_end)
    : allocator(allocator)
    , depth(depth)
    , id_child(id_child)
    , particles_begin(particles_begin)
    , particles_end(particles_end)
    , dim((*particles_begin)->get_position().dim)
    {}
    std::vector<std::vector<std::unique_ptr<NodeI>>> &allocator;
    const unsigned int depth;
    const unsigned int id_child;

    Particle** const particles_begin;
    Particle** const particles_end;

    std::vector<NodeI*> neighbours_list;
    std::vector<NodeI*> interaction_list;

    
    const unsigned int dim;

    Tensor mass_center;
    std::unique_ptr<MultipoleSetI> multipole_set;
    std::unique_ptr<LocalSetI> local_set;

  public:
    Node &get_parent();
    static const std::unique_ptr<MultipoleSetI> &getMultipoleSet(const NodeI &_this); 
    virtual void get_partition(std::vector<std::tuple<Particle*, int, int>> &partitions) const;
    void compute_interaction_list();
    unsigned int get_id() const;
    virtual void compute_multipoles(unsigned int L) = 0;
    virtual void collect_multipoles_to_locals();
    virtual void propagate_locals(LocalSetI *parent_local = nullptr);
  
  protected:
    static const Tensor &getMassCenter(const NodeI &_this); 
    virtual void calculateMC() = 0; //mass center calculator

    friend class Tree;

};

