# Code review feedback

Review scope: C++20 serial and MPI FMM implementation, integration schemes, build configuration, and tests. The source annotations marked `REVIEW` point to the most actionable locations.

## Critical correctness

1. `src/quadratures/runge_kutta_4.cpp:12-36` indexes `accelerations_sum`, `velocity_sum`, and `velocity_k` after `reserve()` only. These vectors have size zero, so any non-empty RK4 step has undefined behaviour. Construct them with `particles.size()` elements (and add a small RK4 regression test).
2. `src/math_utils/tensor.cpp:34-38` uses `reserve()` then `memcpy()` in copy assignment. The destination vector may retain size zero, so the writes are outside its constructed elements; later `operator[]` is also invalid. Replace this body with `data = other.data; dim = data.size();`. The move constructor should likewise move `data` and establish a valid moved-from `dim`.
3. `src/system_updaters/tree_components/tree.hpp:32` dereferences `end` through `&*end`. This is undefined behaviour even for `std::vector` iterators. It also claims generic iterator support while relying on contiguous storage and pointer subtraction. Use a `std::span<const Particle>`/vector interface, or calculate a checked count without dereferencing `end`; explicitly handle empty input.
4. `src/system_updaters/tree_components/parallel_tree.cpp:401-405` gathers computed accelerations and then discards them. Values are packed in unordered-map/Morton traversal order, not original particle order, so they cannot safely be used without particle IDs. Multi-rank states therefore diverge and `ParallelQuadrature` slices by vector index rather than spatial ownership.

## Algorithm and numerical concerns

1. The serial FMM hard-codes order `L = 2` in `src/system_updaters/tree_components/tree.cpp:93-95`; the MPI implementation hard-codes `L = 8`. This makes serial accuracy poor and irreproducible across modes. Expose expansion order, leaf threshold, and softening as validated solver configuration, and add relative-error/error-versus-order tests.
2. `NaiveUpdate` uses the 3D denominator `(r^2 + eps^2) * r` at `src/system_updaters/naive_update.cpp:27-29`, whereas a Plummer-softened gravity law is normally `(r^2 + eps^2)^(3/2)`. FMM leaf P2P terms at `src/system_updaters/tree_components/leaf.cpp:81-95,124-138` are unsoftened. The direct solver is therefore not a consistent reference and coincident distinct particles divide by zero in the FMM path.
3. Tree construction divides by `leaf_size[i]` at `src/system_updaters/tree_components/tree.cpp:49,58`; a system with zero extent along any axis causes division by zero. The depth and Morton-coordinate calculations also use variable-length arrays (`unsigned int n_coor[dim]`), which are a GNU extension and not ISO C++20. Use `std::vector`/fixed dimension types and pad degenerate bounding boxes.
4. The claimed O(N) complexity is plausible only for a bounded leaf occupancy and fixed expansion order. Rebuilding the full tree, allocating many small vectors/maps, and repeatedly allocating expansion objects on every force update add substantial overhead. Profile after correctness fixes; contiguous level storage, preallocated coefficient buffers, and reuse across steps are likely high-impact improvements.

## Organization and build/test quality

1. `CMakeLists.txt:10-20` forces `-O0 -g` for normal builds twice, including `Release`; this prevents meaningful performance measurement. Scope coverage flags to the coverage option and use target-specific compile options. MPI/HDF5 are also required for every build even though serial code/tests do not use MPI; make parallel I/O targets optional or provide a serial configuration.
2. Tests primarily compare single random configurations using absolute tolerances (`tests/test_fmm_vs_naive_2d.cpp`, `tests/test_fmm_vs_naive_3d.cpp`) and emit extensive debug output. They do not cover empty/singleton inputs, repeated/coincident coordinates, degenerate axes, Tensor assignment, RK4, conservation properties, or MPI synchronization. Use deterministic property/regression tests with relative-and-absolute tolerances and retain small diagnostic fixtures only on failure.
3. Parallel ownership, MPI lifecycle, and computation are coupled (`src/system_updaters/parallel_fmm.cpp`, `parallel_tree.cpp`, `parallel_quadrature.hpp`). `ParallelFMM` initializes MPI internally but never finalizes it; application-level code should own MPI initialization/finalization. Separate domain decomposition, communication serialization, and FMM passes so they can be tested independently.

## Validation performed

`cmake -S . -B /tmp/nbody-review-build -DCMAKE_BUILD_TYPE=Release` was attempted. Configuration could not complete because this environment lacks MPI (`Could NOT find MPI_C/MPI_CXX`), so the existing test suite was not runnable here. This is a build-environment limitation, not evidence that the tests pass.
