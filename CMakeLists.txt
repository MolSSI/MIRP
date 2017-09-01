cmake_minimum_required(VERSION 3.2)
cmake_policy(SET CMP0048 NEW) # project VERSION* variables populated from VERSION given below
project(MIRP
        VERSION 0.1
        LANGUAGES C CXX)

enable_testing()

# Standard install locations and cmake helpers
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
set(PN ${PROJECT_NAME})
set(CMAKECONFIG_INSTALL_DIR "share/cmake/mirp")


# Find the dependencies
include(cmake/mirpDependencies.cmake)

# Openmp is optional
set(MIRP_OPENMP CACHE BOOL "Enable OpenMP for some loops")
if(MIRP_OPENMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()


set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS Off)
set(CMAKE_C_EXTENSIONS Off)
add_definitions("-Wall")
add_definitions("-Wextra")
add_definitions("-pedantic")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(mirp)
add_subdirectory(mirp_bin)
add_subdirectory(examples)
add_subdirectory(tests)

##########################################
# Installation and config files
##########################################
# Write out the cmake config file
configure_package_config_file(cmake/mirpConfig.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/mirpConfig.cmake"
                              INSTALL_DESTINATION ${CMAKECONFIG_INSTALL_DIR}
)

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/mirpConfigVersion.cmake
                                 VERSION ${mirp_VERSION}
                                 COMPATIBILITY AnyNewerVersion
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mirpConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/mirpConfigVersion.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR}
)

install(EXPORT "mirpTargets"
        NAMESPACE "mirp::"
        DESTINATION ${CMAKECONFIG_INSTALL_DIR}
)

# Install the dependency-finding script
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/cmake/mirpDependencies.cmake
        DESTINATION ${CMAKECONFIG_INSTALL_DIR}
)

