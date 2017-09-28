cmake_minimum_required(VERSION 3.2)

# Read in the version from the VERSION file
file(STRINGS "VERSION" MIRP_VERSION)

cmake_policy(SET CMP0048 NEW) # project VERSION* variables populated from VERSION given below
project(MIRP
        VERSION ${MIRP_VERSION}
        LANGUAGES C CXX)

enable_testing()

# Standard install locations and cmake helpers
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
set(PN ${PROJECT_NAME})
set(CMAKECONFIG_INSTALL_DIR "share/cmake/mirp")


# add lib64, etc, to search paths (rather than just plain lib)
set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS True)


# Find the dependencies
include(cmake/mirpDependencies.cmake)

# Statically link in as much as possible
set(MIRP_STATIC False CACHE BOOL "Statically link in external libraries as much as possible (libgcc, libstdc++, etc)")

# Openmp is optional
set(MIRP_OPENMP False CACHE BOOL "Enable OpenMP for some loops")
if(MIRP_OPENMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif()


set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS Off)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS Off)


################################
# C COMPILE FLAGS
################################
if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")

    # CMake doesn't support Intel C standard until cmake 3.6
    if("${CMAKE_VERSION}" VERSION_LESS "3.6")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c${CMAKE_C_STANDARD}")
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w3")

    # Globally-disabled warnings
    #   981  : operands are evaluated in unspecified order
    # 11074: Inlining inhibited by limit max-size
    # 11076: To get full report use -qopt-report=4 -qopt-report-phase ipo
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -wd981 -wd11074 -wd11076")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "GNU") 

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -pedantic")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wfloat-equal -Wshadow")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wswitch-default -Wconversion")

elseif("${CMAKE_C_COMPILER_ID}" MATCHES "Clang")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Weverything")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-vla -Wno-padded")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-documentation -Wno-documentation-unknown-command")

endif()


################################
# C++ COMPILE FLAGS
################################
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")

    # CMake doesn't support Intel CXX standard until cmake 3.6
    if("${CMAKE_VERSION}" VERSION_LESS "3.6")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CMAKE_CXX_STANDARD}")
    endif()

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w3")

    # Globally-disabled warnings
    #   981  : operands are evaluated in unspecified order
    #   383: value copied to temporary, reference to temporary used
    # 11074: Inlining inhibited by limit max-size
    # 11076: To get full report use -qopt-report=4 -qopt-report-phase ipo
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd981 -wd383 -wd11074 -wd11076")

    #  2338  : switch statement does not have a default clause
    #list(APPEND CMAKE_CXX_FLAGS "-ww2338")

    if(MIRP_STATIC)
        message(FATAL_ERROR "MIRP_STATIC not available with the Intel compiler")
    endif()

elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfloat-equal -Wshadow")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wswitch-default -Wconversion")

    # We disable sign conversion - there are lots of false positives when
    # using STL containers
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-conversion")

    if(MIRP_STATIC)
        set(CMAKE_EXE_LINKER_FLAGS "-static-libstdc++ -static-libgcc")
    endif()

elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-c++98-compat -Wno-padded")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-documentation -Wno-documentation-unknown-command")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-exit-time-destructors -Wno-global-constructors")

    # We disable sign conversion - there are lots of false positives when
    # using STL containers
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-conversion")

    if(MIRP_STATIC)
        message(FATAL_ERROR "MIRP_STATIC not available with the Clang compiler")
    endif()

endif()


include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(mirp)
add_subdirectory(mirp_bin)
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

# Install the test files
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests
        DESTINATION ${CMAKE_INSTALL_DATADIR}/mirp
)

