cmake_minimum_required(VERSION 3.23)

add_executable(
  ed-test
  test-state.cpp
  test-combinations.cpp
  test-binomial.cpp
  test-hamiltonian.cpp
  test-momentum-space.cpp
)

target_compile_options(
  ed-test
  PUBLIC
  "${ED_CXX_COMPILER_OPTIONS}"
)

target_include_directories(
  ed-test
  PRIVATE
  ${Eigen_INCLUDE_DIRS}
  .
)

target_link_libraries(
    ed-test
    PRIVATE
    gmock
    gmock_main
    gtest
    ed-lib
)

add_test(
    NAME ed-test
    COMMAND ed-test
)
