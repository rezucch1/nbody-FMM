#pragma once
#include <memory>
#include <vector>
#include <tuple>

#include "node.hpp"

class Tree{
  public :
    template <class InputIt>
    Tree(InputIt begin, InputIt end){init_tree(&*begin, &*end);};
    Tensor* get_accelerations(std::vector<Tensor> &acceleration_vector, const Particle* begin, unsigned int size);

  protected:
    std::vector<const Particle*> particle_ordering;
    std::vector<std::vector<std::unique_ptr<NodeI>>> nodes_vector;
    // Node root;

    void init_tree(const Particle* begin, const Particle* end);
  public:
    void print_root_multipoles();
    std::vector<std::tuple<const Particle*, int, int>> get_partition();
    std::vector<std::tuple<NodeI*, NodeI*, unsigned int, char, const Particle*>> get_nodes_interactions();
};