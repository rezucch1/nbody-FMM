#include "system_updaters/tree_components/tree.hpp"
#include <cstdlib>
#include <ctime>

#include "quadratures/symplectic_euler.hpp"
#include "system_updaters/FMM.hpp"

#include <fstream>
// #include <vector>

constexpr double TOTAL_TIME = 100;
constexpr double TIMESTAMP_STEP = 1.0;
std::unique_ptr<QuadratureMethod> quadrature_method(new SyplecticEuler());
std::unique_ptr<SystemUpdateMethod> update_method(new FMM());
constexpr char OUTPUT_FILE[] = "output-FMM-1.0.csv";

constexpr unsigned int RANDOM_SEED = 1234;
constexpr unsigned int NUM_PARTICLES = 15000;
constexpr double DOMAIN_BOUND = 100.0;
constexpr double MAX_MASS = 100.0;
constexpr double MAX_VEL = 100.0;

int main(){
  std::vector<Particle> particles;
  particles.reserve(NUM_PARTICLES);
  srand(RANDOM_SEED);
  for (unsigned int p = 0; p < NUM_PARTICLES; ++p)
    particles.push_back(Particle(rand()*MAX_MASS/RAND_MAX, {rand()*(2*DOMAIN_BOUND)/RAND_MAX - DOMAIN_BOUND, rand()*(2*DOMAIN_BOUND)/RAND_MAX - DOMAIN_BOUND}, {rand()*(2*MAX_VEL)/RAND_MAX - MAX_VEL, rand()*(2*MAX_VEL)/RAND_MAX - MAX_VEL}));

  
  double time = 0.0;
  while (time <= TOTAL_TIME){
    quadrature_method->integrate(update_method.get(), particles, TIMESTAMP_STEP);

    time += TIMESTAMP_STEP;
  }
  // Tree t(particles.begin(), particles.end());
  // auto partitions = t.get_partition();

  // std::ofstream out_file("partitions.csv");

  // out_file << "particle_idx, level, partition_id, x, y" << std::endl;
  // for (const auto &p : partitions)
  //   out_file << std::get<0>(p) - particles.data() << ", "
  //     << std::get<1>(p) << ", "
  //     << std::get<2>(p) << ", "
  //     << std::get<0>(p)->get_position()[0] << ", "
  //     << std::get<0>(p)->get_position()[1] << std::endl;

  // out_file = std::ofstream ("interactions.csv");
  // out_file << "target, source, depth, type, x, y" << std::endl;
  // for (const auto &[target, source, depth, type, particle] : t.get_nodes_interactions()){
  //   out_file << target->get_id() << ", "
  //     << source->get_id() << ", "
  //     << depth << ", "
  //     << type << ", "
  //     << particle->get_position()[0] << ", "
  //     << particle->get_position()[1] << std::endl;
  // }

}



