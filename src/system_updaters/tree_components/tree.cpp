#include "tree.hpp"
#include <tuple>

Tree::Tree(std::vector<Particle> &particles)
: root(nullptr, particles.begin(), particles.end())
{}

