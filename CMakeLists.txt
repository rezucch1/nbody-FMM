cmake_minimum_required(VERSION 3.12.0)
project(nbody LANGUAGES CXX)
cmake_minimum_required()

add_executable(nbody
  src/main.cpp
  src/tensor.cpp
  src/particle.cpp

  src/quadrature_method.cpp
  src/symplectic_euler.cpp
  src/runge_kutta_4.cpp
  
  src/naive_update.cpp
  )
  # src/tree.cpp
  # src/multipole.cpp