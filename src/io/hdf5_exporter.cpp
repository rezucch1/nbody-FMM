/**
 * @file hdf5_exporter.cpp
 * @brief Implementation of HDF5Exporter methods.
 */

#include "hdf5_exporter.hpp"
#include <iostream>
#include <sys/stat.h>

HDF5Exporter::HDF5Exporter(const std::string& dir, const std::string& base)
    : output_dir(dir), base_name(base) {
    mkdir(output_dir.c_str(), 0755);
}

void HDF5Exporter::write_frame(const std::vector<Particle>& particles, int step, double time) {
    std::ostringstream ss;
    ss << base_name << "_" << std::setfill('0') << std::setw(5) << step << ".h5";
    std::string filename = ss.str();
    std::string filepath = output_dir + "/" + filename;

    size_t n = particles.size();
    if (n == 0) return;

    std::vector<double> pos(n * 3, 0.0);
    std::vector<double> vel(n * 3, 0.0);
    std::vector<double> acc(n * 3, 0.0);
    std::vector<double> weight(n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        const auto& p = particles[i];
        pos[i * 3 + 0] = p.get_position()[0];
        pos[i * 3 + 1] = p.get_position()[1];
        pos[i * 3 + 2] = (p.get_position().dim > 2) ? p.get_position()[2] : 0.0;

        vel[i * 3 + 0] = p.get_velocity()[0];
        vel[i * 3 + 1] = p.get_velocity()[1];
        vel[i * 3 + 2] = (p.get_velocity().dim > 2) ? p.get_velocity()[2] : 0.0;

        acc[i * 3 + 0] = p.get_acceleration()[0];
        acc[i * 3 + 1] = p.get_acceleration()[1];
        acc[i * 3 + 2] = (p.get_acceleration().dim > 2) ? p.get_acceleration()[2] : 0.0;

        weight[i] = p.get_weight();
    }

    hid_t file_id = H5Fcreate(filepath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        std::cerr << "Error creating HDF5 file: " << filepath << std::endl;
        return;
    }

    auto write_dataset_2d = [&](const char* name, const std::vector<double>& data) {
        hsize_t dims[2] = {n, 3};
        hid_t space_id = H5Screate_simple(2, dims, NULL);
        hid_t dset_id = H5Dcreate2(file_id, name, H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Dclose(dset_id);
        H5Sclose(space_id);
    };

    auto write_dataset_1d = [&](const char* name, const std::vector<double>& data) {
        hsize_t dims[1] = {n};
        hid_t space_id = H5Screate_simple(1, dims, NULL);
        hid_t dset_id = H5Dcreate2(file_id, name, H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Dclose(dset_id);
        H5Sclose(space_id);
    };

    write_dataset_2d("position", pos);
    write_dataset_2d("velocity", vel);
    write_dataset_2d("acceleration", acc);
    write_dataset_1d("weight", weight);

    H5Fclose(file_id);

    timesteps.push_back(time);
    h5_filenames.push_back(filename);

    write_xdmf_metadata(n);
}

void HDF5Exporter::write_xdmf_metadata(size_t num_particles) {
    std::string xmf_path = output_dir + "/" + base_name + ".xmf";
    std::ofstream xmf(xmf_path);
    if (!xmf.is_open()) return;

    xmf << "<?xml version=\"1.0\" ?>\n";
    xmf << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    xmf << "<Xdmf Version=\"2.0\">\n";
    xmf << " <Domain>\n";
    xmf << "  <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

    for (size_t t = 0; t < timesteps.size(); ++t) {
        xmf << "   <Grid Name=\"Step_" << t << "\" GridType=\"Uniform\">\n";
        xmf << "     <Time Value=\"" << timesteps[t] << "\"/>\n";
        xmf << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << num_particles << "\"/>\n";
        xmf << "     <Geometry GeometryType=\"XYZ\">\n";
        xmf << "       <DataItem Dimensions=\"" << num_particles << " 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
        xmf << "         " << h5_filenames[t] << ":/position\n";
        xmf << "       </DataItem>\n";
        xmf << "     </Geometry>\n";

        xmf << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
        xmf << "       <DataItem Dimensions=\"" << num_particles << " 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
        xmf << "         " << h5_filenames[t] << ":/velocity\n";
        xmf << "       </DataItem>\n";
        xmf << "     </Attribute>\n";

        xmf << "     <Attribute Name=\"Acceleration\" AttributeType=\"Vector\" Center=\"Node\">\n";
        xmf << "       <DataItem Dimensions=\"" << num_particles << " 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
        xmf << "         " << h5_filenames[t] << ":/acceleration\n";
        xmf << "       </DataItem>\n";
        xmf << "     </Attribute>\n";

        xmf << "     <Attribute Name=\"Weight\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        xmf << "       <DataItem Dimensions=\"" << num_particles << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
        xmf << "         " << h5_filenames[t] << ":/weight\n";
        xmf << "       </DataItem>\n";
        xmf << "     </Attribute>\n";

        xmf << "   </Grid>\n";
    }

    xmf << "  </Grid>\n";
    xmf << " </Domain>\n";
    xmf << "</Xdmf>\n";
}
