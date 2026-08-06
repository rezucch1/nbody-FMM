#pragma once
#include <memory>

#include "node.hpp"

class Tree{
  public : Tree(std::vector<Particle> &);
  protected:
    Node root;
};