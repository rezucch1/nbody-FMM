#pragma once
#include <memory>
#include "power_set.hpp"
class NodeI;
#include "node.hpp"

//we need an Interface, with Multipole as a subclass (same with Leaf)
class NodeI{
  protected:
    NodeI(const Node *parent) : parent(parent) {}
    const Node *parent; //pointer to parent multipole (both multipole and leaf have a parent)
    Tensor mass_center;
    std::unique_ptr<MultipoleSetI> multipole_set;

    static const std::unique_ptr<MultipoleSetI> &getMultipoleSet(const NodeI &_this); 
    static const Tensor &getMassCenter(const NodeI &_this); 
    
    virtual void calculateMC() = 0; //mass center calculator

};

