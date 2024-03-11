cmake_minimum_required(VERSION 3.23)

add_executable(
    ed-test
    test-ed.cpp
)

target_include_directories(ed-test PRIVATE .)

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
