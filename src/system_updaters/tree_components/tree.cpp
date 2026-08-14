#include "tree.hpp"
#include <tuple>

#include <iostream>

Tree::Tree(std::vector<Particle> &particles)
: root(nullptr, particles.begin(), particles.end())
{}

void Tree::print_highest_multipole(){
  MultipoleSetI *multipole = NodeI::getMultipoleSet(root).get();

  std::cout << *multipole << std::endl;
}
