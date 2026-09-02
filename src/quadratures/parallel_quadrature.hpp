#pragma once

#include "quadrature_method.hpp"
#include "system_updaters/system_update_method.hpp"
#include <mpi.h>
#include <vector>
#include <type_traits>

/**
 * @file parallel_quadrature.hpp
 * @brief Template wrapper class for ODE integration methods to integrate only a local particle slice per MPI rank.
 */

/**
 * @class NoOpUpdateMethod
 * @brief Internal no-op SystemUpdateMethod passed to inner quadratures when accelerations are pre-computed.
 */
class NoOpUpdateMethod : public SystemUpdateMethod {
  public:
    virtual void update(const std::vector<Particle> &particles) override {
        // Accelerations already calculated and synchronized by ParallelTree
    }
};

/**
 * @class ParallelQuadrature
 * @brief Template wrapper converting any QuadratureMethod into an MPI-parallel slice integrator.
 * @tparam QuadratureType Underlying integration method class (e.g. SymplecticEuler).
 */
template <typename QuadratureType>
class ParallelQuadrature : public QuadratureMethod {
    static_assert(std::is_base_of<QuadratureMethod, QuadratureType>::value,
                  "QuadratureType must inherit from QuadratureMethod");

  public:
    /**
     * @brief Constructs ParallelQuadrature wrapper.
     */
    ParallelQuadrature() : inner_quadrature() {
        int initialized = 0;
        MPI_Initialized(&initialized);
        if (!initialized) {
            rank = 0;
            num_ranks = 1;
        } else {
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
        }
    }

    /**
     * @brief Destructor.
     */
    virtual ~ParallelQuadrature() = default;

    /**
     * @brief Integrates particle positions and velocities over a time step dt on a local particle slice.
     * Delegates acceleration calculations to ParallelTree/Method, then synchronizes updated particle
     * positions and velocities across all MPI ranks using MPI_Allgatherv.
     * @param method Pointer to SystemUpdateMethod solver.
     * @param particles Global vector of particles to advance in time.
     * @param delta_t Time step duration dt.
     */
    virtual void integrate(SystemUpdateMethod *method, std::vector<Particle> &particles, const double &delta_t) override {
        if (particles.empty()) return;

        size_t total_particles = particles.size();
        size_t local_start = (rank * total_particles) / num_ranks;
        size_t local_end = ((rank + 1) * total_particles) / num_ranks;
        size_t local_count = local_end - local_start;

        unsigned int dim = particles[0].get_position().dim;

        // 1. Calculate and synchronize accelerations via ParallelTree/Method
        if (method) {
            method->update(particles);
        }

        // 2. Slice local particles and integrate using inner quadrature with NoOpUpdateMethod
        std::vector<Particle> local_slice(particles.begin() + local_start, particles.begin() + local_end);
        NoOpUpdateMethod noop_solver;
        inner_quadrature.integrate(&noop_solver, local_slice, delta_t);

        // 3. Pack updated local position and velocity tensors into send buffer
        size_t fields_per_particle = 2 * dim; // pos + vel
        std::vector<double> send_buffer;
        send_buffer.reserve(local_count * fields_per_particle);

        for (const auto &p : local_slice) {
            for (unsigned int d = 0; d < dim; ++d) {
                send_buffer.push_back(p.get_position()[d]);
            }
            for (unsigned int d = 0; d < dim; ++d) {
                send_buffer.push_back(p.get_velocity()[d]);
            }
        }

        // 4. Exchange updated particle slice states across all MPI ranks via MPI_Allgatherv
        if (num_ranks > 1) {
            int send_count = (int)send_buffer.size();
            std::vector<int> recv_counts(num_ranks, 0);
            MPI_Allgather(&send_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

            std::vector<int> displs(num_ranks, 0);
            int total_fields = 0;
            for (int r = 0; r < num_ranks; ++r) {
                displs[r] = total_fields;
                total_fields += recv_counts[r];
            }

            std::vector<double> recv_buffer(total_fields, 0.0);
            MPI_Allgatherv(send_buffer.data(), send_count, MPI_DOUBLE,
                           recv_buffer.data(), recv_counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

            // Unpack updated particle states into global particle array
            size_t p_idx = 0;
            size_t i = 0;
            while (i < recv_buffer.size() && p_idx < total_particles) {
                for (unsigned int d = 0; d < dim; ++d) {
                    position(particles[p_idx])[d] = recv_buffer[i++];
                }
                for (unsigned int d = 0; d < dim; ++d) {
                    velocity(particles[p_idx])[d] = recv_buffer[i++];
                }
                p_idx++;
            }
        } else {
            // Single rank: unpack local_slice directly into global particles array
            for (size_t p_idx = local_start; p_idx < local_end; ++p_idx) {
                position(particles[p_idx]) = local_slice[p_idx - local_start].get_position();
                velocity(particles[p_idx]) = local_slice[p_idx - local_start].get_velocity();
            }
        }
    }

  private:
    QuadratureType inner_quadrature;
    int rank;
    int num_ranks;
};
