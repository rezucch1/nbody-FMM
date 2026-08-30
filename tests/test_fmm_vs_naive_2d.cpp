/**
 * @file test_fmm_vs_naive_2d.cpp
 * @brief GoogleTest unit test comparing 2D Naive update method with 2D Fast Multipole Method (FMM).
 */

#include <gtest/gtest.h>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "particle.hpp"
#include "system_updaters/naive_update.hpp"
#include "system_updaters/FMM.hpp"

TEST(FMMvsNaive2DTest, CompareAccelerations2D) {
    constexpr unsigned int NUM_PARTICLES = 20;
    constexpr double DOMAIN_BOUND = 50.0;
    constexpr double PARTICLE_MASS = 1.0e9;
    constexpr unsigned int SEED = 42;

    srand(SEED);

    // 1. Generate N particles in 2D space
    std::vector<Particle> particles_naive;
    particles_naive.reserve(NUM_PARTICLES);

    for (unsigned int i = 0; i < NUM_PARTICLES; ++i) {
        double x = (rand() * 2.0 * DOMAIN_BOUND) / RAND_MAX - DOMAIN_BOUND;
        double y = (rand() * 2.0 * DOMAIN_BOUND) / RAND_MAX - DOMAIN_BOUND;
        double vx = (rand() * 2.0) / RAND_MAX - 1.0;
        double vy = (rand() * 2.0) / RAND_MAX - 1.0;

        particles_naive.push_back(Particle(PARTICLE_MASS, {x, y}, {vx, vy}));
    }

    // Duplicate particle list for FMM comparison
    std::vector<Particle> particles_fmm = particles_naive;

    // 2. Compute accelerations using Naive direct P2P method
    NaiveUpdate naive_solver;
    naive_solver.update(particles_naive);

    // 3. Compute accelerations using Fast Multipole Method (FMM)
    FMM fmm_solver;
    fmm_solver.update(particles_fmm);

    // 4. Compare accelerations particle-by-particle
    std::cout << std::setw(10) << "Particle"
              << std::setw(20) << "Naive Acc Norm"
              << std::setw(20) << "FMM Acc Norm"
              << std::setw(20) << "Abs Error" << std::endl;
    std::cout << std::string(70, '-') << std::endl;

    double max_abs_error = 0.0;

    for (unsigned int i = 0; i < NUM_PARTICLES; ++i) {
        Tensor acc_naive = particles_naive[i].get_acceleration();
        Tensor acc_fmm = particles_fmm[i].get_acceleration();

        double norm_naive = acc_naive.norm();
        double norm_fmm = acc_fmm.norm();

        double abs_error = (acc_naive - acc_fmm).norm();
        if (abs_error > max_abs_error) max_abs_error = abs_error;

        std::cout << std::setw(10) << i
                  << std::setw(20) << std::scientific << norm_naive
                  << std::setw(20) << std::scientific << norm_fmm
                  << std::setw(20) << std::scientific << abs_error << std::endl;
    }

    std::cout << std::string(70, '-') << std::endl;
    std::cout << "[INFO] Maximum acceleration error between 2D Naive and FMM: " << max_abs_error << std::endl;

    // Assert that FMM acceleration closely matches Naive acceleration
    for (unsigned int i = 0; i < NUM_PARTICLES; ++i) {
        Tensor acc_naive = particles_naive[i].get_acceleration();
        Tensor acc_fmm = particles_fmm[i].get_acceleration();

        for (unsigned int d = 0; d < 2; ++d) {
            EXPECT_NEAR(acc_fmm[d], acc_naive[d], 1.0e-3);
        }
    }
}
