# Fast Multipole Method (FMM) for N-Body Problems

[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![HDF5](https://img.shields.io/badge/IO-HDF5%2FXDMF-orange.svg)](https://www.hdfgroup.org/)
[![Apptainer](https://img.shields.io/badge/Container-Singularity%2FApptainer-purple.svg)](https://apptainer.org/)

An advanced, high-performance C++20 implementation of the **Fast Multipole Method (FMM)** for 2D and 3D N-body gravitational and electrostatic particle simulations. This project features exact near-field direct summation, high-order far-field multipole expansions, multiple ODE time-integration schemes, and **HDF5/XDMF** visualization output compatible with ParaView and VisIt.

### Authors & Team Members (Politecnico di Milano)
- **Annea Reyes**
- **Andrei Tomita**
- **Rebecca Zucchi**

**Supervised by**: Prof. Luca Formaggia  
*Advanced Methods for Scientific Computing, September 2026*

> [!NOTE]
> For complete mathematical derivations, spherical harmonic recurrences, and translation theorems, please refer to the technical document:
> **[`docs/FastMultipoleMethodsNbody.pdf`](docs/FastMultipoleMethodsNbody.pdf)**

---

## Table of Contents

- [Overview & Key Features](#overview--key-features)
- [Algorithmic Overview](#algorithmic-overview)
- [Build and Installation](#build-and-installation)
  - [Prerequisites](#prerequisites)
  - [Building Natively with CMake](#building-natively-with-cmake)
  - [Building and Running with Singularity / Apptainer](#building-and-running-with-singularity--apptainer)
- [Running Simulations](#running-simulations)
  - [Running the Main Simulation](#running-the-main-simulation)
  - [Running Unit & Regression Tests](#running-unit--regression-tests)
- [Visualization with ParaView / VisIt](#visualization-with-paraview--visit)
- [References](#references)

---

## Overview & Key Features

Evaluating direct pairwise interactions in an $N$-body system requires $\mathcal{O}(N^2)$ computational complexity, which quickly becomes computationally intractable for large particle numbers. The **Fast Multipole Method (FMM)**, originally developed by Greengard & Rokhlin, reduces this cost to $\mathcal{O}(N)$ or $\mathcal{O}(N \log N)$ by clustering well-separated particles into spatial tree cells and approximating far-field interactions using multipole and local expansion series.

### Key Capabilities

- **2D & 3D FMM Solvers**: Supports 2D Quadtree (complex logarithmic kernel) and 3D Octree (spherical harmonics & solid harmonics) expansion kernels.
- **Morton Z-Order Bit Interleaving**: Fast spatial cell indexing and tree traversal using 64-bit Morton codes (`Morton::encode`, `Morton::decode`).
- **Time Integrators**:
  - **Symplectic Euler** (`SyplecticEuler`): First-order energy-preserving integrator for Hamiltonian systems.
  - **Runge-Kutta 4th Order** (`RungeKutta4`): Classical 4-stage high-accuracy integrator.
- **Direct vs FMM Solvers**:
  - `NaiveUpdate`: Direct $\mathcal{O}(N^2)$ pairwise force solver for baseline comparisons.
  - `FMM`: $\mathcal{O}(N)$ Fast Multipole Method solver.
- **HDF5 / XDMF Exporter** (`HDF5Exporter`): High-throughput time-series output for visualization.

---

## Algorithmic Overview

The FMM algorithm partitions the computational domain into a hierarchical tree (Quadtree in 2D, Octree in 3D) and separates particle interactions into two regimes:

1. **Near-field Interactions**: Nearby particles within adjacent leaf cells are evaluated directly via pairwise summation.
2. **Far-field Interactions**: Well-separated particles are approximated through a series of multipole and local expansions across 4 main passes:
   - **M2M (Multipole-to-Multipole)**: Translates child multipole coefficients up to parent cell centers (Upward Pass).
   - **M2L (Multipole-to-Local)**: Converts far-field multipoles of well-separated source cells in the interaction list into local expansions at target cells (Interaction Phase).
   - **L2L (Local-to-Local)**: Propagates local expansions from parent cells down to child cells (Downward Pass).
   - **Evaluation**: Evaluates local expansions at particle positions to determine far-field forces.

For full mathematical equations, series expansions, and translation matrices, see **[`docs/FastMultipoleMethodsNbody.pdf`](docs/FastMultipoleMethodsNbody.pdf)**.

---

## Build and Installation

### Prerequisites

- **C++ Compiler**: GCC 13+ or Clang 16+ with full C++20 support (`-std=c++20`).
- **Build System**: CMake 3.28+ and Ninja or Make.
- **Data Libraries**: HDF5 (`libhdf5-dev`).
- **Testing Framework**: GoogleTest (fetched automatically via CMake `FetchContent`).

---

### Building Natively with CMake

```bash
# 1. Configure build directory
cmake -B build -DCMAKE_BUILD_TYPE=Release

# 2. Compile executables
cmake --build build -j$(nproc)
```

---

### Building and Running with Singularity / Apptainer

A complete reproducible environment definition is provided in [`Apptainer.def`](Apptainer.def).

```bash
# 1. Build SIF container image
apptainer build --fakeroot build/nbody.sif Apptainer.def

# 2. Build project inside container
apptainer exec build/nbody.sif cmake -B build
apptainer exec build/nbody.sif cmake --build build -j$(nproc)

# 3. Execute unit tests inside container
apptainer exec build/nbody.sif ctest --test-dir build --output-on-failure
```

---

## Running Simulations

### Running the Main Simulation

Runs a 3D FMM simulation with 500 particles and outputs HDF5/XDMF trajectory files:

```bash
./build/nbody
```

---

### Running Unit & Regression Tests

The test suite consists of automated tests covering vector math, direct update correctness, Morton tree construction, and 2D/3D FMM vs Naive accuracy:

```bash
ctest --test-dir build --output-on-failure
```

---

## Visualization with ParaView / VisIt

Simulations generated by `nbody` create an XDMF metadata file and an HDF5 data archive in the `output/` directory:

- `output/nbody_fmm.xmf` (Metadata file referencing HDF5 frames)
- `output/nbody_fmm.h5` (Binary dataset containing particle positions, velocities, and step timestamps)

### Opening Trajectories in ParaView

1. Launch **ParaView**: `paraview`
2. Open File: Select `output/nbody_fmm.xmf`
3. Select Reader: Choose **XDMF Reader**
4. Click **Apply**
5. Add Glyph / Representation: Select **Glyph** $\to$ Glyph Type: **Sphere**
6. Click **Play** to animate particle trajectories over time.
