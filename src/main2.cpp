#include "system_updaters/tree_components/tree.hpp"
#include <cstdlib>
#include <ctime>

#include "quadratures/symplectic_euler.hpp"
#include "system_updaters/naive_update.hpp"
#include "system_updaters/FMM.hpp"
#include "io/hdf5_exporter.hpp"

#include <fstream>
#include <iostream>


#include <complex>
#include <iomanip>

/**
 * @file main2.cpp
 * @brief N-Body simulation driver using FMM and HDF5/XDMF visualization exporter.
 */

constexpr double TOTAL_TIME = 100.0;
constexpr double TIMESTAMP_STEP = 0.1;
constexpr unsigned int RANDOM_SEED = 1234;
constexpr unsigned int NUM_PARTICLES = 100;
constexpr double DOMAIN_BOUND = 100.0;
constexpr double MAX_MASS = 1.0e10;
constexpr double MAX_VEL = 0.0;

constexpr unsigned int DIM = 3;

std::unique_ptr<QuadratureMethod> quadrature_method(new SyplecticEuler());
std::unique_ptr<SystemUpdateMethod> update_method(new FMM());
//std::unique_ptr<SystemUpdateMethod> update_method(new NaiveUpdate());
int main(){
  std::vector<Particle> particles;
  particles.reserve(NUM_PARTICLES);
  srand(RANDOM_SEED);

  for (unsigned int p = 0; p < NUM_PARTICLES; ++p)
{
    if (DIM == 2)
    {
        particles.push_back(
            Particle(
                rand() * MAX_MASS / RAND_MAX,
                {
                    rand() * (2 * DOMAIN_BOUND) / RAND_MAX - DOMAIN_BOUND,
                    rand() * (2 * DOMAIN_BOUND) / RAND_MAX - DOMAIN_BOUND
                },
                {
                    rand() * (2 * MAX_VEL) / RAND_MAX - MAX_VEL,
                    rand() * (2 * MAX_VEL) / RAND_MAX - MAX_VEL
                }
            )
        );
    }
    else if (DIM == 3)
    {
        particles.push_back(
            Particle(
                rand() * MAX_MASS / RAND_MAX,
                {
                    rand() * (2 * DOMAIN_BOUND) / RAND_MAX - DOMAIN_BOUND,
                    rand() * (2 * DOMAIN_BOUND) / RAND_MAX - DOMAIN_BOUND,
                    rand() * (2 * DOMAIN_BOUND) / RAND_MAX - DOMAIN_BOUND
                },
                {
                    rand() * (2 * MAX_VEL) / RAND_MAX - MAX_VEL,
                    rand() * (2 * MAX_VEL) / RAND_MAX - MAX_VEL,
                    rand() * (2 * MAX_VEL) / RAND_MAX - MAX_VEL
                }
            )
        );
    }
}// ============================================================
// TEST get_gradient() in 3D
// ============================================================

/*{
    constexpr unsigned int L = 3;

    LocalSet<3> local(L);

    // Test point
    Tensor d{0.7, -0.4, 0.5};

    // ------------------------------------------------------------
    // Test:
    //
    // L_2^0 = 1 * R_2^0
    //
    // With the chosen normalization:
    //
    // R_2^0 = z^2 - (x^2 + y^2)/2
    //
    // therefore:
    //
    // grad(R_2^0) = (-x, -y, 2z)
    //
    // At d = (0.7, -0.4, 0.5):
    //
    // grad(R_2^0) = (-0.7, 0.4, 1.0)
    // ------------------------------------------------------------


    local.set_elements(
        2,
        0,
        std::complex<double>(1.0, 0.0)
    );

    Tensor gradient = local.get_gradient(d);

    std::cout << "\n========================================\n";
    std::cout << "TEST get_gradient() 3D\n";
    std::cout << "========================================\n";

    std::cout << std::setprecision(15);

    std::cout << "d = ("
              << d[0] << ", "
              << d[1] << ", "
              << d[2] << ")\n";

    std::cout << "Computed gradient:\n";
    std::cout << "("
              << gradient[0] << ", "
              << gradient[1] << ", "
              << gradient[2] << ")\n";

  const double h = 1e-6;

Tensor dx = d;
Tensor dy = d;
Tensor dz = d;

dx[0] += h;
dy[1] += h;
dz[2] += h;

PowerSet<3> p(L, d);
PowerSet<3> px(L, dx);
PowerSet<3> py(L, dy);
PowerSet<3> pz(L, dz);

double gx =
    (px(2, 0).real() - p(2, 0).real()) / h;

double gy =
    (py(2, 0).real() - p(2, 0).real()) / h;

double gz =
    (pz(2, 0).real() - p(2, 0).real()) / h;

std::cout << "Numerical gradient:\n";
std::cout << "("
          << gx << ", "
          << gy << ", "
          << gz << ")\n";
 // R_2^0 = z^2 -(x^2 + y^2) / 2
    // grad(R_2^0) = (-x, -y, 2z)
//
// At d = (0.7, -0.4, 0.5):
// grad(R_2^0) = (-0.7, 0.4, 1.0)

std::cout << "Expected:\n";
std::cout << "(-0.7, 0.4, 1.0)\n";
    std::cout << "========================================\n\n";
}*/
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
