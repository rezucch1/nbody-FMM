/**
 * @file test_naive_update.cpp
 * @brief GoogleTest unit test for NaiveUpdate N-body solver using a 2-particle system.
 */

#include <gtest/gtest.h>
#include <vector>
#include <memory>
#include "particle.hpp"
#include "system_updaters/naive_update.hpp"
#include "quadratures/symplectic_euler.hpp"

TEST(NaiveUpdateTest, TwoParticlesGetCloserUnderGravity) {
    // 1. Setup 2 massive particles at rest along the X axis
    double mass = 1.0e9; // 10^9 kg
    Particle p1(mass, {-10.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    Particle p2(mass, {+10.0, 0.0, 0.0}, {0.0, 0.0, 0.0});

    std::vector<Particle> particles = { p1, p2 };

    double initial_distance = (particles[1].get_position() - particles[0].get_position()).norm();
    EXPECT_DOUBLE_EQ(initial_distance, 20.0);

    // 2. Initialize NaiveUpdate solver and SymplecticEuler integrator
    NaiveUpdate solver;
    SyplecticEuler integrator;

    double dt = 0.5; // Timestep (seconds)
    int num_steps = 20;

    // 3. Integrate over time and verify distance decreases
    double prev_distance = initial_distance;

    for (int step = 1; step <= num_steps; ++step) {
        integrator.integrate(&solver, particles, dt);

        double current_distance = (particles[1].get_position() - particles[0].get_position()).norm();

        // Distance between mutual attractors MUST strictly decrease!
        EXPECT_LT(current_distance, prev_distance);

        // Particle 1 (left) moves right (+x)
        EXPECT_GT(particles[0].get_position()[0], -10.0);

        // Particle 2 (right) moves left (-x)
        EXPECT_LT(particles[1].get_position()[0], +10.0);

        prev_distance = current_distance;
    }

    std::cout << "[INFO] Initial distance: " << initial_distance << " m" << std::endl;
    std::cout << "[INFO] Distance after " << num_steps * dt << " s: " << prev_distance << " m" << std::endl;
}
