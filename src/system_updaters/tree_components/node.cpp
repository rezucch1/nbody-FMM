/**
 * @file node.cpp
 * @brief Implementation of internal Node class methods.
 */

#include "node.hpp"
#include "leaf.hpp"

Node::Node(std::vector<std::unordered_map<uint64_t, std::unique_ptr<NodeI>>> &allocator, unsigned int depth, unsigned int id_child, unsigned int dim)
    : NodeI(allocator, depth, id_child, dim)
{
}

void Node::compute_multipoles(unsigned int L){

    for (auto &c : get_children())
      c.compute_multipoles(L);

    calculateMC();
    if (dim == 2)
      multipole_set = std::make_unique<MultipoleSet<2>>(L);
    else if (dim == 3)
      multipole_set = std::make_unique<MultipoleSet<3>>(L);

    (*multipole_set) = 0;
    for (auto &c : get_children())
      *multipole_set += *getMultipoleSet(c)->weigh_children_with_distance(NodeI::getMassCenter(c) - mass_center);
}

void Node::collect_multipoles_to_locals(){
    NodeI::collect_multipoles_to_locals();

    for (auto &c : get_children())
      c.collect_multipoles_to_locals();
}

void Node::propagate_locals(const LocalSetI *parent_local){
  NodeI::propagate_locals(parent_local);

    for (auto &c : get_children())
      if (local_set) {
        auto child_shift = local_set->distribute_parent_with_distance(getMassCenter(c) - mass_center);
        c.propagate_locals(child_shift.get());
      } else {
        c.propagate_locals(nullptr);
      }
}

const std::vector<NodeI *> &Node::get_neighbours() const
{
  return neighbours_list;
}

NodeI *Node::get_child(unsigned int idx)
{
  auto it = allocator[depth + 1].find((static_cast<uint64_t>(id_child) << dim) + idx);
  return (it != allocator[depth + 1].end()) ? it->second.get() : nullptr;
}

Node::ChildrenContainer Node::get_children(){
  return ChildrenContainer(*this);
}

Node::ChildIterator Node::ChildrenContainer::begin() const {
  auto it = ChildIterator{_this.allocator[_this.depth + 1], _this.id_child << _this.dim, (_this.id_child + 1) << _this.dim};
  it.first();
  return it;
}

Node::ChildIterator Node::ChildrenContainer::end() const {
  return ChildIterator{_this.allocator[_this.depth + 1], (_this.id_child + 1) << _this.dim, (_this.id_child + 1) << _this.dim};
}

void Node::calculateMC()
{
    mass_center = Tensor(dim);
    mass_center *= 0;
    double total_mass = 0;
    for(auto &c : get_children()){
        double child_mass = (*getMultipoleSet(c))(0).real();
        total_mass += child_mass;
        mass_center += child_mass * NodeI::getMassCenter(c);
    }
    mass_center /= total_mass;
    for (int i = 0; i < dim; ++i)
      if (std::isnan(mass_center[i]))
        printf("break\n");

    return;
}
