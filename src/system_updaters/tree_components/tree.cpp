#include "tree.hpp"
#include <tuple>
#include <iterator>

#include <iostream>

void Tree::init_tree(Particle* begin, Particle* end){
  particle_ordering.reserve(end - begin);

  auto p = begin;

  Tensor a = p->get_position();
  Tensor b = p->get_position();

  for (; p < end; ++p){
    particle_ordering.push_back(&*p);
    for (unsigned int d = 0; d < a.dim; ++d){
        if (p->get_position()[d] < a[d])
            a[d] = p->get_position()[d];
        else if (p->get_position()[d] > b[d])
            b[d] = p->get_position()[d];
    }
  }

  auto &root = nodes_vector.emplace_back(1);
  nodes_vector[0][0] = std::make_unique<Node>(nodes_vector, 0, 0, &*particle_ordering.begin(), &*particle_ordering.end(), a, b);

  for (auto &d : nodes_vector)
    for (auto &n : d){
      if (n)
        n->compute_interaction_list();
    }
  
  nodes_vector[0][0]->compute_multipoles(2);
  nodes_vector[0][0]->collect_multipoles_to_locals();
  nodes_vector[0][0]->propagate_locals();
}

void Tree::print_root_multipoles(){
  MultipoleSetI *multipole = NodeI::getMultipoleSet(*nodes_vector[0][0]).get();

  std::cout << *multipole << std::endl;
}

std::vector<std::tuple<Particle *, int, int>> Tree::get_partition()
{
  std::vector<std::tuple<Particle *, int, int>> partitions;
  nodes_vector[0][0]->get_partition(partitions);
  return partitions;
}

std::vector<std::tuple<NodeI *, NodeI *, unsigned int, char, Particle*>> Tree::get_nodes_interactions()
{
  std::vector<std::tuple<NodeI *, NodeI *, unsigned int, char, Particle*>> interaction;
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
