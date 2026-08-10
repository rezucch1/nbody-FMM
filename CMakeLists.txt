cmake_minimum_required(VERSION 3.12.0)
project(nbody LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O0 -g)
include_directories(src)

add_executable(nbody
  src/main2.cpp
)

add_subdirectory(src)

target_link_libraries(nbody PRIVATE
    nbody_core
    math_utils
    quadratures
    system_updaters
    multipole_utils
    tree_components
)