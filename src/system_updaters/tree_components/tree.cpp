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
  unsigned int n_leafs;
  unsigned int depth = 0;
  unsigned int n_leafs_dim;
  Tensor leaf_size;
  std::vector<unsigned int> p_count;

  do{
    n_leafs = (max_n_particles + threshold) / threshold - 1;
    depth += (std::bit_width(n_leafs) + dim + 1) / dim - 1;
    n_leafs = 1 << (depth * dim);
    n_leafs_dim = 1 << depth;
    leaf_size = (b - a) / n_leafs_dim;

    p_count.reserve(n_leafs);
    p_count.assign(n_leafs, 0);

    max_n_particles = 0;
    for (auto p = begin; p < end; ++p){
      Tensor p_pos = p->get_position() - a;
      unsigned int n_coor[dim];
      for (unsigned int i=0; i < dim; ++i){
        n_coor[i] = p_pos[i] / leaf_size[i];
        if (n_coor[i] >= n_leafs_dim) n_coor[i] = n_leafs_dim - 1;
      }
      unsigned int n_idx = Morton::encode(n_coor, dim);
      ++p_count[n_idx];
      if (max_n_particles < p_count[n_idx]) max_n_particles = p_count[n_idx];
    }
  }while (max_n_particles > threshold);

  const Particle** leaf_begin[n_leafs + 1];
  leaf_begin[0] = particle_ordering.data();
  for (unsigned int i = 0; i < n_leafs; ++i){
    leaf_begin[i + 1] = leaf_begin[i] + p_count[i];
    p_count[i] = 0;
  }
  
  for (auto p = begin; p < end; ++p){
    Tensor p_pos = p->get_position() - a;
    unsigned int n_coor[dim];
    for (unsigned int i=0; i < dim; ++i){
      n_coor[i] = p_pos[i] / leaf_size[i];
      if (n_coor[i] >= n_leafs_dim) n_coor[i] = n_leafs_dim - 1;
    }
    unsigned int n_idx = Morton::encode(n_coor, dim);
    leaf_begin[n_idx][p_count[n_idx]] = p;
    ++p_count[n_idx];
  }

  nodes_vector.reserve(depth);
  {
    unsigned int d;
    for (d = 0; d < depth; ++d){
      nodes_vector.emplace_back().reserve(1 << (d * dim));
      unsigned int n_nodes = 0;
      unsigned int i, j;
      i = 0;
      while (i < n_leafs){
        j = i + (1 << ((depth - d) * dim));
        NodeI *node = nullptr;
        if (leaf_begin[j] - leaf_begin[i] > 0)
            node = new Node(nodes_vector, d, n_nodes, leaf_begin[i], leaf_begin[j]);
        nodes_vector[d].emplace_back(node);
        ++n_nodes;
        i = j;
      }
    }

    nodes_vector.emplace_back().reserve(1 << (d * dim));
    for (unsigned int i = 0; i < n_leafs; ++i){
      NodeI *node = nullptr;
      if (leaf_begin[i+1] - leaf_begin[i] > 0)
        node = new Leaf(nodes_vector, d, i, leaf_begin[i], leaf_begin[i+1]);
      nodes_vector[d].emplace_back(node);
    }
  }

  for (auto &d : nodes_vector)
    for (auto &n : d){
      if (n)
        n->compute_interaction_list();
    }
}

void Tree::compute_accelerations(){
  
  nodes_vector[0][0]->compute_multipoles(2);
  nodes_vector[0][0]->collect_multipoles_to_locals();
  nodes_vector[0][0]->propagate_locals();

  for (const auto &l : *(nodes_vector.end() - 1)){
    if (l)
    ((Leaf*) l.get())->compute_acceleration();
  }
}

void Tree::print_root_multipoles(){
  MultipoleSetI *multipole = NodeI::getMultipoleSet(*nodes_vector[0][0]).get();

  std::cout << *multipole << std::endl;
}

std::vector<std::tuple<const Particle *, int, int>> Tree::get_partition()
{
  std::vector<std::tuple<const Particle *, int, int>> partitions;
  partitions.reserve(nodes_vector.size() * particle_ordering.size());

  for (unsigned int d = 0; d < nodes_vector.size(); ++d)
    for (unsigned int n = 0; n < nodes_vector[d].size(); ++n)
      if (nodes_vector[d][n] != nullptr)
        for (auto p = nodes_vector[d][n]->particles_begin; p < nodes_vector[d][n]->particles_end; ++p)
          partitions.emplace_back(*p, d, n);

  return partitions;
}

std::vector<std::tuple<NodeI *, NodeI *, unsigned int, char, const Particle*>> Tree::get_nodes_interactions()
{
  std::vector<std::tuple<NodeI *, NodeI *, unsigned int, char, const Particle*>> interaction;
  for (int d = 0; d < nodes_vector.size(); ++d)
    for (auto &t : nodes_vector[d])if (t){
      for (auto &n : t->interaction_list) if (n)
        for (auto p = n->particles_begin; p < n->particles_end; ++p)
          interaction.emplace_back(t.get(), n, d, 'i', *p);

      for (auto &n : t->neighbours_list) if (n)
        for (auto p = n->particles_begin; p < n->particles_end; ++p)
          interaction.emplace_back(t.get(), n, d, 'n', *p);

      for (auto p = t->particles_begin; p < t->particles_end; ++p)
        interaction.emplace_back(t.get(), t.get(), d, 'x', *p);
    }
    return interaction;
}
