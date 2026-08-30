/**
 * @file test_tree.cpp
 * @brief GoogleTest unit test for FMM Tree spatial construction using 15 particles.
 */

#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "particle.hpp"
#include "system_updaters/tree_components/tree.hpp"
#include "system_updaters/tree_components/leaf.hpp"

class TreeT : public Tree{
    public:
        template<class InputIt>
        TreeT(InputIt begin, InputIt end)
            :Tree(begin, end){}

        unsigned int get_leafs_depth(){
            return nodes_vector.size() - 1;
        }
        
        std::unordered_map<u_int64_t, std::unique_ptr<NodeI>> &get_nodes(unsigned int depth){
            return nodes_vector[depth];
        }

        NodeI &get_node(unsigned int depth, uint64_t idx){
            auto it = nodes_vector[depth].find(idx);
            assert(it != nodes_vector[depth].end());
            return *it->second;
        }
};

struct LeafT : public Leaf{
    public:
        const std::vector<const Particle *> &get_particles(){
            return particles;
        }
};

TEST(TreeTest, ConstructionWith15Particles) {
    constexpr unsigned int NUM_PARTICLES = 15;
    constexpr double DOMAIN_BOUND = 100.0;
    constexpr double PARTICLE_MASS = 10.0;
    constexpr unsigned int SEED = 42;

    srand(SEED);
    std::vector<Particle> particles;
    particles.reserve(NUM_PARTICLES);

    for (unsigned int i = 0; i < NUM_PARTICLES; ++i) {
        double x = (rand() * 2.0 * DOMAIN_BOUND) / RAND_MAX - DOMAIN_BOUND;
        double y = (rand() * 2.0 * DOMAIN_BOUND) / RAND_MAX - DOMAIN_BOUND;
        double z = (rand() * 2.0 * DOMAIN_BOUND) / RAND_MAX - DOMAIN_BOUND;

        particles.push_back(Particle(PARTICLE_MASS, {x, y, z}, {0.0, 0.0, 0.0}));
    }

    TreeT tree(particles.begin(), particles.end());
    LeafT* l;
    const std::vector<const Particle *> *p;

    l = (LeafT*)&tree.get_node(1, 0);
    p = &l->get_particles();
    EXPECT_EQ((*p)[0], &particles[1]);
    EXPECT_EQ((*p)[1], &particles[4]);
    EXPECT_EQ((*p)[2], &particles[7]);
    EXPECT_EQ((*p)[3], &particles[14]);

    l = (LeafT*)&tree.get_node(1, 1);
    p = &l->get_particles();
    EXPECT_EQ((*p)[0], &particles[11]);

    l = (LeafT*)&tree.get_node(1, 3);
    p = &l->get_particles();
    EXPECT_EQ((*p)[0], &particles[2]);
    EXPECT_EQ((*p)[1], &particles[13]);

    l = (LeafT*)&tree.get_node(1, 4);
    p = &l->get_particles();
    EXPECT_EQ((*p)[0], &particles[0]);
    EXPECT_EQ((*p)[1], &particles[3]);

    l = (LeafT*)&tree.get_node(1, 5);
    p = &l->get_particles();
    EXPECT_EQ((*p)[0], &particles[5]);

    l = (LeafT*)&tree.get_node(1, 6);
    p = &l->get_particles();
    EXPECT_EQ((*p)[0], &particles[10]);
    EXPECT_EQ((*p)[1], &particles[12]);

    l = (LeafT*)&tree.get_node(1, 7);
    p = &l->get_particles();
    EXPECT_EQ((*p)[0], &particles[6]);
    EXPECT_EQ((*p)[1], &particles[8]);
    EXPECT_EQ((*p)[2], &particles[9]);

}
