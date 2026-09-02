# Code review feedback

Review scope: C++20 serial and MPI FMM implementation, integration schemes, build configuration, and tests. The source annotations marked `REVIEW` point to the most actionable locations.

## Critical correctness

1. `src/quadratures/runge_kutta_4.cpp:12-36` indexes `accelerations_sum`, `velocity_sum`, and `velocity_k` after `reserve()` only. These vectors have size zero, so any non-empty RK4 step has undefined behaviour. Construct them with `particles.size()` elements (and add a small RK4 regression test), or use emplace_back to add elements to those vectors, not the [] operator. The code may work becouse no bound check is performed (c++23 will change that!)
2. The same problem here! `reserve()` only allocate the buffer but the vector is still logically empty (size=0). `src/math_utils/tensor.cpp:34-38` uses `reserve()` then `memcpy()` in copy assignment. The destination vector may retain size zero, so the writes are outside its constructed elements; later `operator[]` is also invalid. Replace this body with `data = other.data; dim = data.size();`. The move constructor should likewise move `data` and establish a valid moved-from `dim`.
3. `src/system_updaters/tree_components/tree.hpp:32` dereferences `end` through `&*end`. This is undefined behaviour even for `std::vector` iterators. Moreovere you are not using iterators but ponter, which for standard container is discouraged becouse error prone. It also claims generic iterator support while relying on contiguous storage and pointer subtraction.

## Algorithm and numerical concerns

1. The serial FMM hard-codes order `L = 2` in `src/system_updaters/tree_components/tree.cpp:93-95`; the MPI implementation hard-codes `L = 8`. This makes serial accuracy poor and irreproducible across modes. Expose expansion order, leaf threshold, and softening as validated solver configuration, and add relative-error/error-versus-order tests.
2. `NaiveUpdate` uses the 3D denominator `(r^2 + eps^2) * r` at `src/system_updaters/naive_update.cpp:27-29`, whereas a Plummer-softened gravity law is normally `(r^2 + eps^2)^(3/2)`. FMM leaf P2P terms at `src/system_updaters/tree_components/leaf.cpp:81-95,124-138` are unsoftened. The direct solver is therefore not a consistent reference and coincident distinct particles divide by zero in the FMM path.
3. Tree construction divides by `leaf_size[i]` at `src/system_updaters/tree_components/tree.cpp:49,58`; a system with zero extent along any axis causes division by zero. The depth and Morton-coordinate calculations also use variable-length arrays (`unsigned int n_coor[dim]`), which are a GNU extension and not ISO C++20. Use `std::vector`/fixed dimension types and pad degenerate bounding boxes. Did you get the Morton stuff from an external code?


## Organization and build/test quality

1. `CMakeLists.txt:10-20` forces `-O0 -g` for normal builds twice, including `Release`; this prevents meaningful performance measurement. Scope coverage flags to the coverage option and use target-specific compile options. MPI/HDF5 are also required for every build even though serial code/tests do not use MPI; make parallel I/O targets optional or provide a serial configuration.
2. Tests primarily compare single random configurations using absolute tolerances (`tests/test_fmm_vs_naive_2d.cpp`, `tests/test_fmm_vs_naive_3d.cpp`) and emit extensive debug output. They do not cover empty/singleton inputs, repeated/coincident coordinates, degenerate axes, Tensor assignment, RK4... test covering is insufficient.
