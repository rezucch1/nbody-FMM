include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

# Tensor test executable
add_executable(test_tensor test_tensor.cpp)
target_link_libraries(test_tensor PRIVATE math_utils gtest_main)
gtest_discover_tests(test_tensor)

# NaiveUpdate 2-particle test executable
add_executable(test_naive_update test_naive_update.cpp)
target_link_libraries(test_naive_update PRIVATE 
    nbody_core 
    system_updaters 
    tree_components 
    multipole_utils 
    quadratures 
    math_utils 
    gtest_main
)
gtest_discover_tests(test_naive_update)

# 3D Naive vs FMM comparison test executable
add_executable(test_fmm_vs_naive_3d test_fmm_vs_naive_3d.cpp)

target_link_libraries(test_fmm_vs_naive_3d PRIVATE 
    nbody_core 
    system_updaters 
    tree_components 
    multipole_utils 
    quadratures 
    math_utils 
    gtest_main
)

gtest_discover_tests(test_fmm_vs_naive_3d)

# Tree construction 15-particle test executable
add_executable(test_tree test_tree.cpp)
target_link_libraries(test_tree PRIVATE 
    nbody_core 
    system_updaters 
    tree_components 
    multipole_utils 
    quadratures 
    math_utils 
    gtest_main
)
gtest_discover_tests(test_tree)

# 2D Naive vs FMM comparison test executable
add_executable(test_fmm_vs_naive_2d test_fmm_vs_naive_2d.cpp)
target_link_libraries(test_fmm_vs_naive_2d PRIVATE 
    nbody_core 
    system_updaters 
    tree_components 
    multipole_utils 
    quadratures 
    math_utils 
    gtest_main
)
gtest_discover_tests(test_fmm_vs_naive_2d)
