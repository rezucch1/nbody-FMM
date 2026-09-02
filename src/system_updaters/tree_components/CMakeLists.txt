add_library(tree_components OBJECT
    leaf.cpp
    node_i.cpp
    node.cpp
    tree.cpp
    parallel_tree.cpp
    ghost_node.cpp
    ghost_leaf.cpp
    morton_code.cpp
)

target_include_directories(tree_components PUBLIC
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/math_utils
    ${CMAKE_SOURCE_DIR}/src/system_updaters/multipole_utils
    ${CMAKE_SOURCE_DIR}/src/system_updaters/tree_components
)

target_link_libraries(tree_components PUBLIC MPI::MPI_CXX)
