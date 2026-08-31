/**
 * @file ghost_node.cpp
 * @brief Implementation of GhostNode class.
 */

#include "ghost_node.hpp"

GhostNode::GhostNode(std::vector<std::unordered_map<uint64_t, std::unique_ptr<NodeI>>> &allocator, unsigned int depth, uint64_t id_child, unsigned int dim)
: GhostNodeI(allocator, depth, id_child, dim)
{}
