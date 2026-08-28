#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <hdf5.h>
#include "particle.hpp"

/**
 * @file hdf5_exporter.hpp
 * @brief Exporter for writing N-Body particle trajectories to HDF5 binary files and XDMF time series for ParaView.
 */

/**
 * @class HDF5Exporter
 * @brief Class managing per-timestamp HDF5 binary dataset creation and XDMF metadata indexing.
 */
class HDF5Exporter {
  private:
    std::string output_dir;
    std::string base_name;
    std::vector<double> timesteps;
    std::vector<std::string> h5_filenames;

  public:
    /**
     * @brief Constructs HDF5Exporter with destination directory and base filename.
     * @param dir Output directory path (e.g. "output").
     * @param base File base prefix (e.g. "particles").
     */
    HDF5Exporter(const std::string& dir = "output", const std::string& base = "particles");

    /**
     * @brief Exports a single timestep frame of particles to an HDF5 binary file and updates XDMF index.
     * @param particles Const reference to particle list.
     * @param step Integer timestep index.
     * @param time Simulation time value.
     */
    void write_frame(const std::vector<Particle>& particles, int step, double time);

  private:
    /**
     * @brief Writes master .xmf XML index file connecting all HDF5 timestamps for ParaView.
     * @param num_particles Number of particles in dataset.
     */
    void write_xdmf_metadata(size_t num_particles);
};
