add_library(io
    hdf5_exporter.cpp
)
target_include_directories(io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${HDF5_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(io PUBLIC ${HDF5_LIBRARIES} MPI::MPI_CXX)
