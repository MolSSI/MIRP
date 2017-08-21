cmake_minimum_required(VERSION 3.0)
project(Arbint LANGUAGES C CXX)
enable_testing()


macro(mirp_find_library libname)
    find_library(MIRP_${libname}_LIBRARY ${libname})

    if(NOT MIRP_${libname}_LIBRARY)
        message(FATAL_ERROR "Library ${libname} not found!")
    endif()

    get_filename_component(MIRP_${libname}_DIR ${MIRP_${libname}_LIBRARY} DIRECTORY)
    get_filename_component(MIRP_${libname}_DIR ${MIRP_${libname}_DIR} DIRECTORY)
    message(STATUS "Found ${libname} library: ${MIRP_${libname}_LIBRARY}")
    message(STATUS "Found ${libname} path: ${MIRP_${libname}_DIR}")
    include_directories(${MIRP_${libname}_DIR}/include)
endmacro()

add_definitions("-Wall")
add_definitions("-Wextra")
add_definitions("-pedantic")

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS Off)
set(CMAKE_C_EXTENSIONS Off)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
mirp_find_library(arb)
mirp_find_library(flint)
mirp_find_library(mpfr)
mirp_find_library(gmp)

add_subdirectory(mirp)
add_subdirectory(mirp_bin)
add_subdirectory(examples)
add_subdirectory(tests)
