cmake_minimum_required(VERSION 3.23)

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Spectra REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")
find_package(OpenBLAS REQUIRED)


add_executable(
  ed
  main.cpp
)

target_compile_options(
  ed
  PUBLIC
  "${ED_CXX_COMPILER_OPTIONS}"
)

set_target_properties(ed PROPERTIES RUNTIME_OUTPUT_DIRECTORY ..)
target_link_libraries(ed PUBLIC ed-lib)

add_library(
    ed-lib
    Timing.cpp
)

target_include_directories(
  ed-lib
  PUBLIC
  /usr/include/eigen3
  /usr/include/Spectra
  ${OpenBLAS_INCLUDE_DIR}
  .
)

target_link_libraries(
    ed-lib
    PUBLIC
    OpenMP::OpenMP_CXX
    ${OpenBLAS_LIB}
    gfortran
)

