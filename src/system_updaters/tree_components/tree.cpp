#include "tree.hpp"
#include <tuple>
#include <iterator>

#include <iostream>

Node Tree::init_root(Tree* _this, Particle* begin, Particle* end){
  _this->particle_ordering.reserve(end - begin);

  auto p = begin;

  Tensor a = p->get_position();
  Tensor b = p->get_position();

  for (; p < end; ++p){
    _this->particle_ordering.push_back(&*p);
    for (unsigned int d = 0; d < a.dim; ++d){
        if (p->get_position()[d] < a[d])
            a[d] = p->get_position()[d];
        else if (p->get_position()[d] > b[d])
            b[d] = p->get_position()[d];
    }
  }

  return Node(nullptr, &*_this->particle_ordering.begin(), &*_this->particle_ordering.end(), a, b);
}

void Tree::print_root_multipoles(){
  MultipoleSetI *multipole = NodeI::getMultipoleSet(root).get();

  std::cout << *multipole << std::endl;
}
