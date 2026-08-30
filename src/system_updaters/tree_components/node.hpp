#pragma once

#include <memory>
#include <list>
#include <unordered_map>
#include <cstdint>
#include "particle.hpp"
#include "math_utils/tensor.hpp"
#include <cassert>

class Node;
#include "node_i.hpp"

/**
 * @file node.hpp
 * @brief Internal non-leaf Node class in the FMM spatial tree representation.
 */

/**
 * @class Node
 * @brief Represents an internal tree node managing child nodes and aggregate multipole expansions.
 */
class Node : public NodeI {
  protected:

  struct ChildrenContainer;

  /**
 * @brief Custom iterator stepping over child nodes in sparse unordered_map allocator.
 */
  struct ChildIterator {
    std::unordered_map<uint64_t, std::unique_ptr<NodeI>> &map;
    uint64_t offset;
    const uint64_t max_offset;
    NodeI *node = nullptr;

    bool operator!=(const ChildIterator& other) const { return offset != other.offset; }
    bool operator<(const ChildIterator& other) const { return offset < other.offset; }
    
    void first(){
      std::unordered_map<uint64_t, std::unique_ptr<NodeI>>::iterator it;
      while(offset < max_offset){
        it = map.find(offset);
        if (it != map.end()){
          node = it->second.get();
          return;
        }
        ++offset;
      }
      node = nullptr;
    }

    ChildIterator& operator++() {
      ++offset;
      first();
      return *this;
    }
    // ChildIterator operator+(size_t n) const { return ChildIterator{map, base_id, offset + n}; }
    
    NodeI &operator*() const {
      assert(node != nullptr);
      return *node;
    }

    NodeI &operator->() const {
      return operator*();
    }
  };

    struct ChildrenContainer{
      Node &_this;
      ChildIterator begin() const;
      ChildIterator end() const;
    };
    

  public:

    /**
     * @brief Constructs an internal Node and initializes child nodes.
     * @param allocator Reference to 2D nodes grid.
     * @param depth Tree level depth.
     * @param id_child Morton child ID.
     * @param dim Spatial dimension.
     */
    Node(std::vector<std::unordered_map<uint64_t, std::unique_ptr<NodeI>>> &allocator, unsigned int depth, unsigned int id_child, unsigned int dim);

    /**
     * @brief Computes multipole expansions recursively upward from children (M2M phase).
     * @param L Maximum degree of expansion.
     */
    virtual void compute_multipoles(unsigned int L) override;

    /**
     * @brief Collects multipoles from interaction list into local expansion set (M2L phase).
     */
    virtual void collect_multipoles_to_locals() override;

    /**
     * @brief Propagates parent local expansion down to children (L2L phase).
     * @param parent_local Translated local expansion from parent node.
     */
    virtual void propagate_locals(const LocalSetI *parent_local = nullptr) override;

    /**
     * @brief Gets list of neighbor nodes.
     * @return Reference to neighbor node pointers vector.
     */
    const std::vector<NodeI *> &get_neighbours() const;

  protected:
    /**
     * @brief Gets child node pointer by index.
     * @param idx Child index (0 to \f$ 2^d - 1 \f$).
     * @return Pointer to child node.
     */
    NodeI* get_child(unsigned int idx);

    ChildrenContainer get_children();

    /**
     * @brief Calculates center of mass from child node centers of mass.
     */
    virtual void calculateMC() override;
};