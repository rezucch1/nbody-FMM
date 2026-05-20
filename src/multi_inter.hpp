#pragma once
#include <memory>
#include "power_set.hpp"
class Multi_inter;
#include "node.hpp"

//we need an Interface, with Multipole as a subclass (same with Leaf)
class Multi_inter{
  protected:
    Multi_inter(const Node *parent) : parent(parent) {}
    const Node *parent; //pointer to parent multipole (both multipole and leaf have a parent)
    Tensor mass_center;
    std::unique_ptr<PowerSetI> multipole_set;
    
    virtual void calculateMC()=0; //mass center calculator

};

