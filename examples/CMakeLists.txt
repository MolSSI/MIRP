cmake_minimum_required(VERSION 3.2)
enable_testing()

# Read in the version from the VERSION file
get_filename_component(PARENT_DIR ${CMAKE_CURRENT_LIST_DIR} DIRECTORY)
file(STRINGS "${PARENT_DIR}/VERSION" MIRP_VERSION)

# add lib64, etc, to search paths (rather than just plain lib)
set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS True)

# Common C/C++ flags
include(${PARENT_DIR}/cmake/mirpFlags.cmake)

#############################
# Find MIRP
#############################
find_package(mirp ${MIRP_VERSION} REQUIRED)

# Example 1 - simple eri
add_executable(simple_eri simple_eri.c)
target_link_libraries(simple_eri mirp::mirp)
add_test(NAME simple_eri COMMAND simple_eri)


# Example 1 - reading a reference file in C++
add_executable(read_ref read_ref.cpp)
add_test(NAME read_ref COMMAND read_ref ${CMAKE_CURRENT_LIST_DIR}/../tests/eri_water_sto-3g.ref)
