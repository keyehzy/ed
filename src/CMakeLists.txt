cmake_minimum_required(VERSION 3.23)

add_executable(
    ed
    main.cpp
)

set_target_properties(ed PROPERTIES RUNTIME_OUTPUT_DIRECTORY ..)
target_link_libraries(ed PRIVATE ed-lib)

add_library(
    ed-lib
    foo.cpp
)
target_include_directories(ed-lib PUBLIC .)

