#include "system_updaters/tree_components/tree.hpp"
#include <cstdlib>
#include <ctime>

#include <fstream>
// #include <vector>

constexpr unsigned int RANDOM_SEED = 1234;
constexpr unsigned int NUM_PARTICLES = 150;
constexpr double DOMAIN_BOUND = 100.0;
constexpr double MAX_MASS = 100.0;
constexpr double MAX_VEL = 100.0;

int main(){
  std::vector<Particle> particles;
  particles.reserve(NUM_PARTICLES);
  srand(RANDOM_SEED);
  for (unsigned int p = 0; p < NUM_PARTICLES; ++p)
    particles.push_back(Particle(rand()*MAX_MASS/RAND_MAX, {rand()*(2*DOMAIN_BOUND)/RAND_MAX - DOMAIN_BOUND, rand()*(2*DOMAIN_BOUND)/RAND_MAX - DOMAIN_BOUND}, {rand()*(2*MAX_VEL)/RAND_MAX - MAX_VEL, rand()*(2*MAX_VEL)/RAND_MAX - MAX_VEL}));
  Tree t(particles.begin(), particles.end());
  auto partitions = t.get_partition();

  std::ofstream out_file("partitions.csv");

  out_file << "particle_idx, level, partition_id, x, y" << std::endl;
  for (const auto &p : partitions)
    out_file << std::get<0>(p) - particles.data() << ", "
      << std::get<1>(p) << ", "
      << std::get<2>(p) << ", "
      << std::get<0>(p)->get_position()[0] << ", "
      << std::get<0>(p)->get_position()[1] << std::endl;

  out_file = std::ofstream ("interactions.csv");
  out_file << "target, source, depth, type, x, y" << std::endl;
  for (const auto &[target, source, depth, type, particle] : t.get_nodes_interactions()){
    out_file << target->get_id() << ", "
      << source->get_id() << ", "
      << depth << ", "
      << type << ", "
      << particle->get_position()[0] << ", "
      << particle->get_position()[1] << std::endl;
  }

}



