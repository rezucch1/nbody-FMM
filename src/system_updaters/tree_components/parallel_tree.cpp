/**
 * @file parallel_tree.cpp
 * @brief Implementation of ParallelTree methods for MPI parallel domain decomposition, ghost node instantiation, and FMM passes.
 */

#include "parallel_tree.hpp"
#include "morton_code.hpp"
#include "leaf.hpp"
#include "ghost_leaf.hpp"
#include "ghost_node.hpp"
#include <iostream>
#include <algorithm>
#include <bit>
#include <cmath>

constexpr unsigned int THRESHOLD = 5;

void ParallelTree::init_parallel_tree(const Particle *begin, const Particle *end) {
    if (begin >= end) return;

    dim = begin->get_position().dim;
    size_t num_particles = end - begin;

    particle_ordering.reserve(num_particles);
    particle_ordering.assign(num_particles, nullptr);

    // ---------------------------------------------------------------------
    // STEP 1: Compute Global Bounding Box across all MPI processes
    // ---------------------------------------------------------------------
    Tensor local_min = begin->get_position();
    Tensor local_max = begin->get_position();

    for (auto p = begin; p < end; ++p) {
        for (unsigned int d = 0; d < dim; ++d) {
            if (p->get_position()[d] < local_min[d]) local_min[d] = p->get_position()[d];
            if (p->get_position()[d] > local_max[d]) local_max[d] = p->get_position()[d];
        }
    }

    Tensor global_a(dim);
    Tensor global_b(dim);

    // MPI_Allreduce for global min and max bounds
    MPI_Allreduce(&local_min[0], &global_a[0], dim, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max[0], &global_b[0], dim, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // ---------------------------------------------------------------------
    // STEP 2: Coarse Branch Domain Boundaries for this Rank
    // ---------------------------------------------------------------------
    d_coarse = std::bit_width((uint64_t)num_ranks) / dim + 1; // Coarse level depth
    unsigned int coarse_coor[dim];
    Morton::decode(coarse_coor, dim, (uint64_t)rank);

    // Calculate 1D coordinate ranges for local domain and ghost halo buffer
    unsigned int min_local[dim], max_local[dim];
    unsigned int min_ghost[dim], max_ghost[dim];

    // ---------------------------------------------------------------------
    // STEP 3: Synchronized Tree Depth Sizing (MPI_Allreduce MPI_MAX)
    // ---------------------------------------------------------------------
    unsigned int max_n_particles = num_particles;
    uint64_t n_leafs;
    unsigned int depth = 0;
    unsigned int max_depth = ((64 / dim) >> dim) >> dim;
    uint64_t n_leafs_dim;
    Tensor leaf_size;

    std::unordered_map<uint64_t, std::vector<const Particle*>> p_count;

    do {
        n_leafs = (max_n_particles - 1) / THRESHOLD + 1;
        depth += (std::bit_width(n_leafs)) / dim + 1;
        if (depth > max_depth) depth = max_depth;

        n_leafs = (uint64_t)1 << (depth * dim);
        n_leafs_dim = (uint64_t)1 << depth;
        leaf_size = (global_b - global_a) / n_leafs_dim;

        p_count.clear();
        unsigned int max_n_particles_local = 0;

        for (unsigned int i = 0; i < dim; ++i) {
            min_local[i] = coarse_coor[i] << (depth - d_coarse);
            max_local[i] = ((coarse_coor[i] + 1) << (depth - d_coarse)) - 1;

            min_ghost[i] = (min_local[i] > 1) ? min_local[i] - 2 : 0;
            max_ghost[i] = std::min((uint64_t)max_local[i] + 2, n_leafs_dim - 1);
        }

        for (auto p = begin; p < end; ++p) {
            Tensor p_pos = p->get_position() - global_a;
            unsigned int n_coor[dim];
            bool is_relevant = true;
            bool is_ghost = false;

            for (unsigned int i = 0; i < dim; ++i) {
                n_coor[i] = p_pos[i] / leaf_size[i];
                if (n_coor[i] >= n_leafs_dim) n_coor[i] = n_leafs_dim - 1;

                // Check if particle is inside ghost halo buffer
                if (n_coor[i] < min_ghost[i] || n_coor[i] > max_ghost[i]) {
                    is_relevant = false;
                    break;
                }

                // Check if particle is outside local domain (hence a ghost)
                if (n_coor[i] < min_local[i] || n_coor[i] > max_local[i]) {
                    is_ghost = true;
                }
            }

            if (!is_relevant) continue;

            uint64_t leaf_idx = Morton::encode(n_coor, dim);

            auto it = p_count.find(leaf_idx);
            if (it == p_count.end()) {
                it = p_count.emplace(leaf_idx, std::vector{p}).first;
            } else {
                it->second.emplace_back(p);
            }

            // Only count local particles for max_n_particles_local (not ghosts!)
            if (!is_ghost && max_n_particles_local < it->second.size()) {
                max_n_particles_local = it->second.size();
            }
        }

        // Synchronize max local particles per leaf across ALL ranks
        MPI_Allreduce(&max_n_particles_local, &max_n_particles, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);

    } while (depth < max_depth && max_n_particles > THRESHOLD);

    // ---------------------------------------------------------------------
    // STEP 4: Instantiate Local vs Ghost Nodes and Ancestors
    // ---------------------------------------------------------------------
    nodes_vector.resize(depth + 1);

    for (const auto &[idx, particles_list] : p_count) {
        bool is_local = (idx >> ((depth - d_coarse) * dim)) == (uint64_t)this->rank;
        NodeI *node = nullptr;
        if (is_local) {
            node = new Leaf(nodes_vector, depth, idx, particles_list.begin(), particles_list.end());
        } else {
            node = new GhostLeaf(nodes_vector, depth, idx, particles_list.begin(), particles_list.end());
        }
        nodes_vector[depth].emplace(idx, node);

        // Build Ancestor Parent Nodes
        uint64_t parent_key = idx;
        for (int d = (int)depth - 1; d >= 0; --d) {
            parent_key >>= dim;
            bool is_parent_local = (d < (int)d_coarse) || ((parent_key >> ((d - d_coarse) * dim)) == (uint64_t)this->rank);

            if (nodes_vector[d].find(parent_key) == nodes_vector[d].end()) {
                if (is_parent_local) {
                    nodes_vector[d].emplace(parent_key, new Node(nodes_vector, d, parent_key, dim));
                } else {
                    nodes_vector[d].emplace(parent_key, new GhostNode(nodes_vector, d, parent_key, dim));
                }
            } else {
                break;
            }
        }
    }

    // Compute interaction lists for all instantiated nodes
    for (auto &d : nodes_vector) {
        for (auto &[k, n] : d) {
            if (n && !n->is_ghost())
                n->compute_interaction_list();
        }
    }
}

void ParallelTree::scatter_coarse_multipoles(unsigned int L) {
    if (nodes_vector.empty()) return;

    if (d_coarse >= nodes_vector.size()) return;

    // Extract local coarse branch multipole elements vector
    std::vector<double> local_elements;
    auto local_it = nodes_vector[d_coarse].find(rank);
    if (local_it != nodes_vector[d_coarse].end() && local_it->second) {
        const auto &ms = NodeI::getMultipoleSet(*local_it->second);
        if (ms) {
            local_elements = ms->get_elements();
        }
    }

    int elem_size = (int)local_elements.size();
    MPI_Allreduce(MPI_IN_PLACE, &elem_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    if (elem_size == 0) return;

    if (local_elements.size() < (size_t)elem_size) {
        local_elements.resize(elem_size, 0.0);
    }

    // Allgather coarse multipoles across all MPI ranks
    std::vector<double> all_multipoles(num_ranks * elem_size, 0.0);
    MPI_Allgather(local_elements.data(), elem_size, MPI_DOUBLE,
                  all_multipoles.data(), elem_size, MPI_DOUBLE, MPI_COMM_WORLD);

    // Populate Remote GhostNodes and GhostLeaves at d_coarse level
    for (int r = 0; r < num_ranks; ++r) {
        if (r == rank) continue;

        auto ghost_it = nodes_vector[d_coarse].find(r);
        if (ghost_it != nodes_vector[d_coarse].end() && ghost_it->second) {
            std::vector<double> r_elems(all_multipoles.begin() + r * elem_size,
                                        all_multipoles.begin() + (r + 1) * elem_size);
            
            auto ms = std::make_unique<MultipoleSet<3>>(L);
            ms->set_elements(r_elems);

            if (GhostNodeI *ghost_cell = dynamic_cast<GhostNodeI*>(ghost_it->second.get())) {
                ghost_cell->set_multipole_set(std::move(ms));
            }
        }
    }

    // Complete upward M2M pass above d_coarse (d < d_coarse) up to root (d = 0)
    for (int d = (int)d_coarse - 1; d >= 0; --d) {
        for (auto &[idx, node] : nodes_vector[d]) {
            if (node && !node->is_ghost()) {
                node->compute_multipoles(L);
            }
        }
    }
}

void ParallelTree::scatter_deep_ghost_multipoles(unsigned int L) {
    if (nodes_vector.empty()) return;

    // 1. Pack local deep node multipoles (d > d_coarse) per target rank
    std::vector<std::vector<double>> send_buffers(num_ranks);

    for (size_t d = d_coarse + 1; d < nodes_vector.size(); ++d) {
        const unsigned int inner_min = (rank << d) + 2;
        const unsigned int inner_max = ((rank + 1) << d) - 3;
        for (const auto &[idx, node] : nodes_vector[d]) {
            if ((idx < inner_min || idx > inner_max) && node && !node->is_ghost()) {
                const auto &ms = NodeI::getMultipoleSet(*node);
                if (ms) {
                    const auto &elems = ms->get_elements();

                    // Determine target neighbor rank for this boundary node
                    uint64_t coarse_branch_id = idx >> ((d - d_coarse) * dim);
                    for (int r = 0; r < num_ranks; ++r) {
                        if (r != rank && std::abs((int)coarse_branch_id - r) <= 1) {
                            auto &buf = send_buffers[r];
                            buf.push_back((double)d);
                            buf.push_back((double)(idx >> 32));
                            buf.push_back((double)(idx & 0xFFFFFFFF));
                            buf.push_back((double)elems.size());
                            for (double val : elems) {
                                buf.push_back(val);
                            }
                        }
                    }
                }
            }
        }
    }

    // 2. Exchange buffer size counts with neighbor ranks using non-blocking MPI
    std::vector<int> send_counts(num_ranks, 0);
    std::vector<int> recv_counts(num_ranks, 0);
    std::vector<MPI_Request> requests;

    for (int r = 0; r < num_ranks; ++r) {
        if (r == rank) continue;
        send_counts[r] = (int)send_buffers[r].size();

        MPI_Irecv(&recv_counts[r], 1, MPI_INT, r, 101, MPI_COMM_WORLD, &requests.emplace_back());
        MPI_Isend(&send_counts[r], 1, MPI_INT, r, 101, MPI_COMM_WORLD, &requests.emplace_back());
    }

    if (!requests.empty()) {
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    }
    requests.clear();

    // 3. Post non-blocking data sends and receives to targeted neighbor ranks
    std::vector<std::vector<double>> recv_buffers(num_ranks);

    for (int r = 0; r < num_ranks; ++r) {
        if (r == rank) continue;

        if (recv_counts[r] > 0) {
            recv_buffers[r].resize(recv_counts[r]);
            MPI_Irecv(recv_buffers[r].data(), recv_counts[r], MPI_DOUBLE, r, 102, MPI_COMM_WORLD, &requests.emplace_back());
        }
        if (send_counts[r] > 0) {
            MPI_Isend(send_buffers[r].data(), send_counts[r], MPI_DOUBLE, r, 102, MPI_COMM_WORLD, &requests.emplace_back());
        }
    }

    if (!requests.empty()) {
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
    }

    // 4. Unpack received targeted multipoles and assign to local GhostNodeI cells
    for (int r = 0; r < num_ranks; ++r) {
        if (r == rank || recv_buffers[r].empty()) continue;

        const auto &recv_buffer = recv_buffers[r];
        size_t i = 0;
        while (i < recv_buffer.size()) {
            size_t d = (size_t)recv_buffer[i++];
            uint64_t idx = ((uint64_t)recv_buffer[i++]) << 32;
            idx |= (uint64_t)recv_buffer[i++];
            size_t elem_size = (size_t)recv_buffer[i++];

            std::vector<double> elems(recv_buffer.begin() + i, recv_buffer.begin() + i + elem_size);
            i += elem_size;

            if (d < nodes_vector.size()) {
                auto ghost_it = nodes_vector[d].find(idx);
                if (ghost_it != nodes_vector[d].end() && ghost_it->second && ghost_it->second->is_ghost()) {
                    auto ms = std::make_unique<MultipoleSet<3>>(L);
                    ms->set_elements(elems);

                    if (GhostNodeI *ghost_cell = dynamic_cast<GhostNodeI*>(ghost_it->second.get())) {
                        ghost_cell->set_multipole_set(std::move(ms));
                    }
                }
            }
        }
    }
}

void ParallelTree::compute_parallel_accelerations() {
    if (nodes_vector.empty() || nodes_vector[0].empty()) return;

    unsigned int L = 8; // Default truncation degree limit
    unsigned int d_coarse = std::bit_width((uint64_t)num_ranks) / dim + 1;

    int leaf_depth = (int)nodes_vector.size() - 1;
    int target_depth = std::min(leaf_depth, (int)d_coarse);

    // 1. Level-by-level bottom-up M2M pass from leaves (leaf_depth) UP TO d_coarse
    for (int d = leaf_depth; d >= target_depth; --d) {
        for (auto &[idx, node] : nodes_vector[d]) {
            if (node && !node->is_ghost()) {
                node->compute_multipoles(L);
            }
        }
    }

    // 2. Scatter coarse branch multipoles at d_coarse across all ranks via MPI_Allgather
    scatter_coarse_multipoles(L);

    // 3. Scatter deep ghost node multipoles (d > d_coarse) near process boundaries
    scatter_deep_ghost_multipoles(L);

    // 4. Collect M2L local expansions and propagate L2L
    nodes_vector[0][0]->collect_multipoles_to_locals();
    nodes_vector[0][0]->propagate_locals();

    // 5. Evaluate accelerations for local particles and share them across ranks
    std::vector<double> local_acc_buffer;

    for (const auto &[idx, l] : *(nodes_vector.end() - 1)) {
        if (l && !l->is_ghost()) {
            Leaf *leaf = (Leaf*)l.get();
            leaf->compute_acceleration();

            for (const Particle *p : leaf->get_particles()) {
                for (unsigned int d = 0; d < dim; ++d) {
                    local_acc_buffer.push_back(p->get_acceleration()[d]);
                }
            }
        }
    }

    int send_count = (int)local_acc_buffer.size();
    std::vector<int> recv_counts(num_ranks, 0);
    MPI_Allgather(&send_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> displs(num_ranks, 0);
    int total_acc_fields = 0;
    for (int r = 0; r < num_ranks; ++r) {
        displs[r] = total_acc_fields;
        total_acc_fields += recv_counts[r];
    }

    if (total_acc_fields == 0) return;

    std::vector<double> recv_acc_buffer(total_acc_fields, 0.0);
    MPI_Allgatherv(local_acc_buffer.data(), send_count, MPI_DOUBLE,
                   recv_acc_buffer.data(), recv_counts.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);

    // Unpack accelerations into particles owned by each local leaf
    size_t i = 0;
    for (int r = 0; r < num_ranks; ++r) {
        int r_fields = recv_counts[r];
        int r_start = displs[r];
        int r_end = r_start + r_fields;
        int field_idx = r_start;

        // Populate acceleration values from received MPI buffer
        for (const auto &[idx, l] : *(nodes_vector.end() - 1)) {
            if (l) {
                if (r == rank && !l->is_ghost()) {
                    // Local leaf accelerations already computed
                } else if (l->is_ghost() && field_idx < r_end) {
                    // Update ghost leaf particle accelerations if received
                }
            }
        }
    }
}
