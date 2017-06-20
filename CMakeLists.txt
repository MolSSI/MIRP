cmake_minimum_required(VERSION 3.0)
project(Arbint LANGUAGES C CXX)

add_definitions("-Wall")
add_definitions("-Wextra")
add_definitions("-pedantic")

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_C_EXTENSIONS ON)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_subdirectory(mirp)
add_subdirectory(mirppy)
add_subdirectory(test)
