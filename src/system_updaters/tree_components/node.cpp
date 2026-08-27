/**
 * @file node.cpp
 * @brief Implementation of internal Node class methods.
 */

#include "node.hpp"
#include "leaf.hpp"
constexpr unsigned int threshold = 5;

Node::Node(std::vector<std::vector<std::unique_ptr<NodeI>>> &allocator, unsigned int depth, unsigned int id_child, const Particle** particles_begin, const Particle** particles_end)
    : NodeI(allocator, depth, id_child ,particles_begin, particles_end)
{
}

void Node::compute_multipoles(unsigned int L){

    for (auto c = get_children_begin(); c < get_children_end(); ++c) if (*c)
        (*c)->compute_multipoles(L);

    calculateMC();

    multipole_set = std::make_unique<MultipoleSet<2>>(L);
    (*multipole_set) = 0;
    for (auto c = get_children_begin(); c < get_children_end(); ++c) if (*c)
        *multipole_set += *getMultipoleSet(**c)->weigh_children_with_distance(NodeI::getMassCenter(**c) - mass_center);
}

void Node::collect_multipoles_to_locals(){
    NodeI::collect_multipoles_to_locals();

    for (auto c = get_children_begin(); c < get_children_end(); ++c) if (*c)
        (*c)->collect_multipoles_to_locals();
}

void Node::propagate_locals(LocalSetI *parent_local){
    NodeI::propagate_locals(parent_local);

    for (auto c = get_children_begin(); c < get_children_end(); ++c) if (*c)
        if (local_set)
          (*c)->propagate_locals(local_set->distribute_parent_with_distance(getMassCenter(**c) - mass_center));
        else 
          (*c)->propagate_locals(nullptr);
}

const std::vector<NodeI *> &Node::get_neighbours() const
{
  return neighbours_list;
}

void Node::get_partition(std::vector<std::tuple<const Particle *, int, int>> &partitions) const{
    NodeI::get_partition(partitions);

    for (auto c = get_children_begin(); c < get_children_end(); ++c) if (*c)
        (*c)->get_partition(partitions);
}

NodeI *Node::get_child(unsigned int idx)
{
  return allocator[depth + 1][id_child * 1<<dim + idx].get();
}

std::vector<std::unique_ptr<NodeI>>::iterator Node::get_children_begin(){
  return allocator[depth + 1].begin() + (id_child << dim);
}

std::vector<std::unique_ptr<NodeI>>::iterator Node::get_children_end(){
  return allocator[depth + 1].begin() + ((id_child+1) << dim);
}

std::vector<std::unique_ptr<NodeI>>::const_iterator Node::get_children_begin() const
{
  return allocator[depth + 1].cbegin() + (id_child << dim);
}

std::vector<std::unique_ptr<NodeI>>::const_iterator Node::get_children_end() const
{
  return allocator[depth + 1].cbegin() + ((id_child+1) << dim);
}

void Node::calculateMC()
{
    mass_center = Tensor(dim);
    mass_center *= 0;
    double total_mass = 0;
    for(auto c = get_children_begin(); c < get_children_end(); ++c) if (*c){
        double child_mass = (*getMultipoleSet(**c))(0).real();
        total_mass += child_mass;
        mass_center = child_mass * NodeI::getMassCenter(**c);
    }
    mass_center /= total_mass;

    return;
}
