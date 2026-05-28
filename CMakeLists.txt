cmake_minimum_required(VERSION 3.12.0)
project(nbody LANGUAGES CXX)
cmake_minimum_required()

# add_executable(nbody
#   src/main.cpp
#   src/tensor.cpp
#   src/particle.cpp

#   src/quadrature_method.cpp
#   src/symplectic_euler.cpp
#   src/runge_kutta_4.cpp
  
#   src/naive_update.cpp
#   )
  # src/tree.cpp
  # src/multipole.cpp

add_executable(nbody
  src/main2.cpp

  src/tensor.cpp
  src/particle.cpp

  src/tree.cpp
  src/node_i.cpp
  src/node.cpp
  src/leaf.cpp

  src/base_multipole_set_i.cpp
  src/multipole_set.cpp
  src/multipole_set2D.cpp
  src/power_set.cpp
  src/power_set2D.cpp

  src/math_utils.hpp
)