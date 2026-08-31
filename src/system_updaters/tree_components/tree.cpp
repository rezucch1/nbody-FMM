/**
 * @file tree.cpp
 * @brief Implementation of Tree spatial partitioning, Morton index sorting, and FMM pass orchestration.
 */

#include "tree.hpp"
#include <tuple>
#include <iterator>
#include <bit>
#include <iostream>

#include "morton_code.hpp"
#include "leaf.hpp"

void Tree::init_tree(const Particle *begin, const Particle *end)
{
  particle_ordering.reserve(end - begin);
  particle_ordering.assign(end - begin, nullptr);

  Tensor a = begin->get_position();
  Tensor b = begin->get_position();

  for (auto p = begin; p < end; ++p){
    for (unsigned int d = 0; d < a.dim; ++d){
      if (p->get_position()[d] < a[d])
          a[d] = p->get_position()[d];
      else if (p->get_position()[d] > b[d])
          b[d] = p->get_position()[d];
    }
  }

  unsigned int threshold = 5;
  unsigned int dim = begin->get_position().dim;

  unsigned int max_n_particles = end - begin;
  uint64_t n_leafs;
  unsigned int depth = 0;
  unsigned int max_depth = ((64 / dim) >> dim) << dim;
  uint64_t n_leafs_dim;
  Tensor leaf_size;
  std::unordered_map<uint64_t, std::vector<const Particle*>> p_count;

  do{
    n_leafs = (max_n_particles - 1) / threshold + 1;
    depth += (std::bit_width(n_leafs)) / dim + 1;
    if (depth > max_depth) depth = max_depth;
    n_leafs = (uint64_t)1 << (depth * dim);
    n_leafs_dim = (uint64_t)1 << depth;
    leaf_size = (b - a) / n_leafs_dim;

    p_count.clear();

    max_n_particles = 0;
    for (auto p = begin; p < end; ++p){
      Tensor p_pos = p->get_position() - a;
      unsigned int n_coor[dim];
      for (unsigned int i=0; i < dim; ++i){
        n_coor[i] = p_pos[i] / leaf_size[i];
        if (n_coor[i] >= n_leafs_dim) n_coor[i] = n_leafs_dim - 1;
      }
      uint64_t n_idx = Morton::encode(n_coor, dim);
      auto it = p_count.find(n_idx);
      if (it == p_count.end())
        it = p_count.emplace(n_idx, std::vector{p}).first;
      else
        it->second.emplace_back(p);
      if (max_n_particles < it->second.size()) max_n_particles = it->second.size();
    }
  }while (depth < max_depth && max_n_particles > threshold);

  nodes_vector.resize(depth + 1);
  for (const auto &[idx, particles] : p_count){
    NodeI *node = new Leaf(nodes_vector, depth, idx, particles.begin(), particles.end());
    nodes_vector[depth].emplace(idx, node);

    uint64_t parent_key = idx;
    for (int d = (int)depth - 1; d >= 0; --d) {
      parent_key >>= dim;
      if (nodes_vector[d].find(parent_key) == nodes_vector[d].end()) {
        nodes_vector[d].emplace(parent_key, new Node(nodes_vector, d, parent_key, begin->get_position().dim));
      } else 
        break;
    }
  }

  for (auto &d : nodes_vector)
    for (auto &[k, n] : d){
      if (n)
        n->compute_interaction_list();
    }
}

void Tree::compute_accelerations(){
  unsigned int L = 8;
  
  for (int d = (int)nodes_vector.size() - 1; d >= 0; --d) {
    for (auto &[idx, node] : nodes_vector[d]) {
      if (node) {
        node->compute_multipoles(L);
      }
    }
  }

  nodes_vector[0][0]->collect_multipoles_to_locals();
  nodes_vector[0][0]->propagate_locals();

  for (const auto &[idx, l] : *(nodes_vector.end() - 1)){
    if (l)
      ((Leaf*) l.get())->compute_acceleration();
  }
}

void Tree::print_root_multipoles(){
  MultipoleSetI *multipole = NodeI::getMultipoleSet(*nodes_vector[0][0]).get();

  std::cout << *multipole << std::endl;
}

// std::vector<std::tuple<const Particle *, int, int>> Tree::get_partition()
// {
//   std::vector<std::tuple<const Particle *, int, int>> partitions;
//   partitions.reserve(nodes_vector.size() * particle_ordering.size());

//   for (unsigned int d = 0; d < nodes_vector.size(); ++d)
//     for (unsigned int n = 0; n < nodes_vector[d].size(); ++n)
//       if (nodes_vector[d][n] != nullptr)
//         for (auto p = nodes_vector[d][n]->particles_begin; p < nodes_vector[d][n]->particles_end; ++p)
//           partitions.emplace_back(*p, d, n);

//   return partitions;
// }

// std::vector<std::tuple<NodeI *, NodeI *, unsigned int, char, const Particle*>> Tree::get_nodes_interactions()
// {
//   std::vector<std::tuple<NodeI *, NodeI *, unsigned int, char, const Particle*>> interaction;
//   for (int d = 0; d < nodes_vector.size(); ++d)
//     for (auto &t : nodes_vector[d])if (t){
//       for (auto &n : t->interaction_list) if (n)
//         for (auto p = n->particles_begin; p < n->particles_end; ++p)
//           interaction.emplace_back(t.get(), n, d, 'i', *p);

//       for (auto &n : t->neighbours_list) if (n)
//         for (auto p = n->particles_begin; p < n->particles_end; ++p)
//           interaction.emplace_back(t.get(), n, d, 'n', *p);

//       for (auto p = t->particles_begin; p < t->particles_end; ++p)
//         interaction.emplace_back(t.get(), t.get(), d, 'x', *p);
//     }
//     return interaction;
// }
