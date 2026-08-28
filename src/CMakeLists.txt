add_subdirectory(math_utils)
add_subdirectory(quadratures)
add_subdirectory(system_updaters)
add_subdirectory(io)

add_library(nbody_core OBJECT
    particle.cpp
)
