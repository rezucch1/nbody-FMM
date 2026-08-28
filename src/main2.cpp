#include "system_updaters/tree_components/tree.hpp"
#include <cstdlib>
#include <ctime>

#include "quadratures/symplectic_euler.hpp"
#include "system_updaters/naive_update.hpp"
#include "system_updaters/FMM.hpp"
#include "io/hdf5_exporter.hpp"

#include <fstream>
#include <iostream>

/**
 * @file main2.cpp
 * @brief N-Body simulation driver using FMM and HDF5/XDMF visualization exporter.
 */

constexpr double TOTAL_TIME = 100.0;
constexpr double TIMESTAMP_STEP = 1.0;
constexpr unsigned int RANDOM_SEED = 1234;
constexpr unsigned int NUM_PARTICLES = 15000;
constexpr double DOMAIN_BOUND = 100.0;
constexpr double MAX_MASS = 100.0;
constexpr double MAX_VEL = 10.0;

std::unique_ptr<QuadratureMethod> quadrature_method(new SyplecticEuler());
std::unique_ptr<SystemUpdateMethod> update_method(new FMM());

int main(){
  std::vector<Particle> particles;
  particles.reserve(NUM_PARTICLES);
  srand(RANDOM_SEED);

  for (unsigned int p = 0; p < NUM_PARTICLES; ++p)
    particles.push_back(Particle(rand()*MAX_MASS/RAND_MAX, {rand()*(2*DOMAIN_BOUND)/RAND_MAX - DOMAIN_BOUND, rand()*(2*DOMAIN_BOUND)/RAND_MAX - DOMAIN_BOUND, rand()*(2*DOMAIN_BOUND)/RAND_MAX - DOMAIN_BOUND}, {rand()*(2*MAX_VEL)/RAND_MAX - MAX_VEL, rand()*(2*MAX_VEL)/RAND_MAX - MAX_VEL, rand()*(2*MAX_VEL)/RAND_MAX - MAX_VEL}));

  HDF5Exporter hdf5_exporter("output", "nbody_fmm");

  double time = 0.0;
  unsigned int step = 0;

  std::cout << "Starting FMM simulation (N=" << NUM_PARTICLES << ") with HDF5/XDMF export..." << std::endl;

  while (time <= TOTAL_TIME){
    hdf5_exporter.write_frame(particles, step, time);
    quadrature_method->integrate(update_method.get(), particles, TIMESTAMP_STEP);
    time += TIMESTAMP_STEP;
    ++step;
  }

  std::cout << "Simulation complete. Trajectory exported to output/nbody_fmm.xmf" << std::endl;

  return 0;
}
