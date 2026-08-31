#include "system_updaters/tree_components/tree.hpp"
#include <cstdlib>
#include <ctime>

#include "quadratures/symplectic_euler.hpp"
#include "system_updaters/naive_update.hpp"
#include "system_updaters/FMM.hpp"
#include "io/hdf5_exporter.hpp"

#include <fstream>
#include <iostream>
#include <cmath>
#include <memory>
#include <numbers>
#include <utility>
#include <vector>

/**
 * @file main2.cpp
 * @brief N-Body simulation driver using FMM and HDF5/XDMF visualization exporter.
 */

constexpr double TOTAL_TIME = 1000.0;
constexpr double TIMESTAMP_STEP = 1;
constexpr unsigned int RANDOM_SEED = 1234;
constexpr unsigned int NUM_PARTICLES = 500;
constexpr double DOMAIN_BOUND = 100.0;
constexpr double MAX_MASS = 1.0e10;
constexpr double MAX_VEL = 0.0;

constexpr unsigned int DIM = 2;

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
}
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

/*

constexpr double TOTAL_TIME = 500.0;
constexpr double TIMESTAMP_STEP = 1.0; // Ridotto per stabilità numerica dell'integratore
constexpr unsigned int RANDOM_SEED = 1234;
constexpr unsigned int NUM_PARTICLES = 500;
constexpr double SPHERE_RADIUS = 100.0;
constexpr double TOTAL_MASS = 1.0e12;   // Massa totale del sistema
constexpr double SOFTENING_EPS = 1.0;    // Previene singolarità 1/r per r->0

std::unique_ptr<QuadratureMethod> quadrature_method(new SyplecticEuler());
std::unique_ptr<SystemUpdateMethod> update_method(new FMM());

// Generatore di numeri casuali uniformi in [0, 1]
double rand_unif() {
    return static_cast<double>(rand()) / RAND_MAX;
}

// Generatore di numeri casuali secondo distribuzione Gaussiana N(0, 1) tramite Box-Muller
double rand_normal() {
    double u1 = rand_unif();
    double u2 = rand_unif();
    while (u1 <= 1e-15) u1 = rand_unif(); // Evita log(0)
    return std::sqrt(-2.0 * std::log(u1)) *
           std::cos(2.0 * std::numbers::pi * u2);
}

int main() {
    std::vector<Particle> particles;
    particles.reserve(NUM_PARTICLES);
    std::vector<Tensor> positions;
    positions.reserve(NUM_PARTICLES);
    std::vector<Tensor> velocities;
    velocities.reserve(NUM_PARTICLES);
    srand(RANDOM_SEED);

    const double particle_mass = TOTAL_MASS / NUM_PARTICLES;

    // 1. GENERAZIONE POSIZIONI IN UN DOMINIO SFERICO E VELOCITÀ CASUALI
    for (unsigned int p = 0; p < NUM_PARTICLES; ++p) {
        // Distribuizione sferica uniforme (Metodo del raggio cubico)
        double r = SPHERE_RADIUS * std::cbrt(rand_unif());
        double theta = std::acos(1.0 - 2.0 * rand_unif()); // Angolo polare
        double phi = 2.0 * std::numbers::pi * rand_unif();  // Angolo azimutale

        double x = r * std::sin(theta) * std::cos(phi);
        double y = r * std::sin(theta) * std::sin(phi);
        double z = r * std::cos(theta);

        // Stima iniziale della dispersione di velocità sigma
        double sigma = std::sqrt(GRAV_CONST * TOTAL_MASS /
                                 (6.0 * SPHERE_RADIUS));

        double vx = sigma * rand_normal();
        double vy = sigma * rand_normal();
        double vz = sigma * rand_normal();

        positions.emplace_back(Tensor{x, y, z});
        velocities.emplace_back(Tensor{vx, vy, vz});
    }

    // 2. CORREZIONE MOTO DEL CENTRO DI MASSA (V_CM = 0)
    double v_cm_x = 0.0, v_cm_y = 0.0, v_cm_z = 0.0;
    for (const auto& velocity : velocities) {
        v_cm_x += velocity[0];
        v_cm_y += velocity[1];
        v_cm_z += velocity[2];
    }
    v_cm_x /= NUM_PARTICLES;
    v_cm_y /= NUM_PARTICLES;
    v_cm_z /= NUM_PARTICLES;

    for (auto& velocity : velocities) {
        velocity[0] -= v_cm_x;
        velocity[1] -= v_cm_y;
        velocity[2] -= v_cm_z;
    }

    // 3. CALCOLO ENERGIA POTENZIALE (U) E CINETICA (K) PER IL VIRIALE
    double U = 0.0;
    double K = 0.0;

    for (size_t i = 0; i < positions.size(); ++i) {
        // Cinetica
        K += 0.5 * particle_mass * velocities[i].squared_norm();

        // Potenziale tra coppie
        for (size_t j = i + 1; j < positions.size(); ++j) {
            const Tensor distance = positions[i] - positions[j];
            const double dist = std::sqrt(distance.squared_norm() +
                                          SOFTENING_EPS * SOFTENING_EPS);

            U -= GRAV_CONST * particle_mass * particle_mass / dist;
        }
    }

    // 4. RISCALAMENTO VIRIALE DELLE VELOCITÀ (Impone 2K + U = 0)
    double alpha = std::sqrt(std::abs(U) / (2.0 * K));
    for (auto& velocity : velocities) {
        velocity *= alpha;
    }

    for (size_t i = 0; i < positions.size(); ++i)
        particles.emplace_back(particle_mass, std::move(positions[i]),
                               std::move(velocities[i]));

    std::cout << "Sistema Virializzato correttamente. Fattore di scala alpha = " << alpha << std::endl;

    // 5. CICLO DI SIMULAZIONE
    HDF5Exporter hdf5_exporter("output", "nbody_fmm");

    double time = 0.0;
    unsigned int step = 0;

    std::cout << "Starting FMM simulation (N=" << NUM_PARTICLES << ") in stable sphere..." << std::endl;

    while (time <= TOTAL_TIME) {
        hdf5_exporter.write_frame(particles, step, time);
        quadrature_method->integrate(update_method.get(), particles, TIMESTAMP_STEP);

        time += TIMESTAMP_STEP;
        ++step;
    }

    std::cout << "Simulation complete. Trajectory exported to output/nbody_fmm.xmf" << std::endl;

    return 0;
}

*/