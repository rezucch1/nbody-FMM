cmake_minimum_required(VERSION 3.12.0)
project(nbody LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

option(ENABLE_COVERAGE "Enable GCC/Clang code coverage reporting" OFF)
if(ENABLE_COVERAGE)
    add_compile_options(--coverage -O0 -g)
    add_link_options(--coverage)
else()
    add_compile_options(-O0 -g)
endif()

find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)

add_compile_options(-O0 -g)
include_directories(src ${HDF5_INCLUDE_DIRS} ${MPI_CXX_INCLUDE_PATH})

add_executable(nbody
  src/main.cpp
)

add_subdirectory(src)
add_subdirectory(tests)

target_link_libraries(nbody PRIVATE
    nbody_core
    io
    math_utils
    quadratures
    system_updaters
    multipole_utils
    tree_components
    ${HDF5_LIBRARIES}
    MPI::MPI_CXX
)