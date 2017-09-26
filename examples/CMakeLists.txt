cmake_minimum_required(VERSION 3.2)
enable_testing()

# Read in the version from the VERSION file
get_filename_component(PARENT_DIR ${CMAKE_CURRENT_LIST_DIR} DIRECTORY)
file(STRINGS "${PARENT_DIR}/VERSION" MIRP_VERSION)

# Compile flags
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS Off)
set(CMAKE_C_EXTENSIONS Off)

#############################
# Find MIRP
#############################
find_package(mirp ${MIRP_VERSION} REQUIRED)

# Example 1 - simple eri
add_executable(simple_eri simple_eri.c)
target_link_libraries(simple_eri mirp::mirp)
add_test(NAME simple_eri COMMAND simple_eri)
