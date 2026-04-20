#include "particle.hpp"
#include "tensor.hpp"

#include "symplectic_euler.hpp"
#include "runge_kutta_4.hpp"
#include "naive_update.hpp"

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <memory>


constexpr unsigned int N_PARTICLES = 4;
constexpr double RADIUS = 1.0;
constexpr double MASS = 1.0;
constexpr double TOTAL_TIME = 10000000;
constexpr double TIMESTAMP_STEP = 1.0;
std::unique_ptr<QuadratureMethod> quadrature_method(new RungeKutta4());
std::unique_ptr<SystemUpdateMethod> update_method(new NaiveUpdate());
constexpr char OUTPUT_FILE[] = "output-RK4-1.0.csv";
constexpr unsigned int MAX_POINTS_OUTPUT = 10000;


int main(){

  double v0 = 0.0;
  for (unsigned int i = 1; i < N_PARTICLES; ++i)
    v0 += 1 / std::sin(M_PI * ((double) i / N_PARTICLES));
  
  v0 = 0.5 * sqrt(GRAV_CONST * MASS * v0 / RADIUS);

  std::vector<Particle> particles;
  particles.reserve(N_PARTICLES);

  for (unsigned int i = 0; i < N_PARTICLES; ++i){
    double alpha =  M_PI * ((double) 2 * i / N_PARTICLES);
    particles.push_back(Particle(MASS, RADIUS * Tensor{std::cos(alpha), std::sin(alpha)}, v0 * Tensor{-std::sin(alpha), std::cos(alpha)}));
  }

  std::ofstream csv_stream_file(OUTPUT_FILE);
  csv_stream_file << "time";
  for (unsigned int i = 0; i < N_PARTICLES; ++i)
    csv_stream_file << ", x" << i << ", y" << i;
  csv_stream_file << std::endl;

  csv_stream_file << "0";
  for (unsigned int i = 0; i < N_PARTICLES; ++i){
    Tensor pos = particles[i].get_position();
    csv_stream_file << ", " << pos[0] << ", " << pos[1];
  }
  csv_stream_file << std::endl;

  double time = 0.0;
  
  unsigned int prev_point_bucket_id = 0;
  unsigned int point_bucket_id = 0;

  while (time <= TOTAL_TIME){
    quadrature_method->integrate(update_method.get(), particles, TIMESTAMP_STEP);

    // Write the points of the current timestamp on output file only if them are from a different
    // bucket respect the previusly timestamp
    prev_point_bucket_id = point_bucket_id;
    point_bucket_id = (unsigned int)((time / TOTAL_TIME) * (MAX_POINTS_OUTPUT / N_PARTICLES));

    if (prev_point_bucket_id != point_bucket_id){
      csv_stream_file << time;
      for (unsigned int i = 0; i < N_PARTICLES; ++i){
        Tensor pos = particles[i].get_position();
        csv_stream_file << ", " << pos[0] << ", " << pos[1];
      }
      csv_stream_file << std::endl;
    }

    time += TIMESTAMP_STEP;
  }

  std::cout << std::endl;

  // std::cout << "Particle ID, Mass, Possition, Velocity, Distance from center" << std::endl;
  // for (unsigned int i = 0; i < N_PARTICLES; i++){
  //   Particle p = particles[i];
  //   std::cout << i << ", " << p.get_mass() << ", " << p.get_position() << ", " << p.get_velocity() << ", " << p.get_position().norm() << std::endl;
  // }

  return 0;
}