add_subdirectory(multipole_utils)
add_subdirectory(tree_components)

add_library(system_updaters OBJECT
    naive_update.cpp
    FMM.cpp
)
