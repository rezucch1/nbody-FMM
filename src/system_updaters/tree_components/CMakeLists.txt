add_library(tree_components OBJECT
    leaf.cpp
    node_i.cpp
    node.cpp
    tree.cpp
    morton_code.cpp
)

target_include_directories(tree_components PUBLIC
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/math_utils
    ${CMAKE_SOURCE_DIR}/src/system_updaters/multipole_utils
    ${CMAKE_SOURCE_DIR}/src/system_updaters/tree_components
)
