#pragma once
#include <memory>
#include <vector>
#include <tuple>

#include "node.hpp"

class Tree{
  public :
    template <class InputIt>
    Tree(InputIt begin, InputIt end){init_tree(&*begin, &*end);};

  protected:
    std::vector<Particle*> particle_ordering;
    std::vector<std::vector<std::unique_ptr<NodeI>>> nodes_vector;
    // Node root;

    void init_tree(Particle* begin, Particle* end);
  public:
    void print_root_multipoles();
    std::vector<std::tuple<Particle*, int, int>> get_partition();
    std::vector<std::tuple<NodeI*, NodeI*, unsigned int, char, Particle*>> get_nodes_interactions();
};