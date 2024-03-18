cmake_minimum_required(VERSION 3.23)

project(ed CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CXX_FLAGS "-Wall -Wextra -g")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(ED_CXX_COMPILER_OPTIONS "" CACHE STRING "")
mark_as_advanced(ED_CXX_COMPILER_OPTIONS)

add_subdirectory(src)
add_subdirectory(third-party)

if (BUILD_TESTING)
  include(CTest)
  add_subdirectory(test)
endif()
