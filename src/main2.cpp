#include "tree.hpp"
#include "cstdlib"
#include "ctime"
// #include <vector>

constexpr unsigned int NUM_PARTICLES = 100;
constexpr double DOMAIN_BOUND = 100.0;
constexpr double MAX_MASS = 100.0;
constexpr double MAX_VEL = 100.0;

int main(){
  std::vector<Particle> particles;
  particles.reserve(NUM_PARTICLES);
  srand(time(NULL));
  for (unsigned int p = 0; p < NUM_PARTICLES; ++p)
    particles.push_back(Particle(rand()*MAX_MASS/RAND_MAX, {rand()*(2*DOMAIN_BOUND)/RAND_MAX - DOMAIN_BOUND, rand()*(2*DOMAIN_BOUND)/RAND_MAX - DOMAIN_BOUND}, {rand()*(2*MAX_VEL)/RAND_MAX - MAX_VEL, rand()*(2*MAX_VEL)/RAND_MAX - MAX_VEL}));
  Tree t(particles);
}



