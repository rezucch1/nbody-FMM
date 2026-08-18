#pragma once
#include <memory>
#include <vector>
#include <tuple>

#include "node.hpp"

class Tree{
  public :
    template <class InputIt>
    Tree(InputIt begin, InputIt end) : root(init_root(this, &*begin, &*end)){}

  protected:
    std::vector<Particle*> particle_ordering; 
    Node root;

    static Node init_root(Tree* _this, Particle* begin, Particle* end);
  public:
    void print_root_multipoles();
    std::vector<std::tuple<Particle*, int, int>> get_partition();
};