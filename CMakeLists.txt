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

# CMake doesn't support Intel CXX standard until cmake 3.6
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
  if("${CMAKE_VERSION}" VERSION_LESS "3.6")
      add_compile_options(-std=c++${CMAKE_CXX_STANDARD})
  endif()
endif()


if("${CMAKE_C_COMPILER_ID}" MATCHES "Intel")
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
endif()


if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w3")

    # Globally-disabled warnings
    #   981  : operands are evaluated in unspecified order
    #   383: value copied to temporary, reference to temporary used
    # 11074: Inlining inhibited by limit max-size
    # 11076: To get full report use -qopt-report=4 -qopt-report-phase ipo
    set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -wd981 -wd383 -wd11074 -wd11076")

    #  2338  : switch statement does not have a default clause
    #list(APPEND CMAKE_CXX_FLAGS "-ww2338")

    # some useless remarks
    #  1418  : external function definition with no prior declaration
    #  1419  : external declaration in primary source file 
    #   383  : value copied to temporary, reference to temporary used (generally ok in c++)
    #   981  : operands are evaluated in unspecified order (generally ok in c++)
    #   869  : parameter "XXX" was never referenced.
    #          This is due to a bug in Intel before Intel 16 that causes this
    #          remark to be printed for variadic templates
    # 11074  : Inlining inhibited by limit max-size (etc)
    # 11076  : To get full report use -qopt-report=4 -qopt-report-phase ipo
    #   444  : destructor for base class .. not virtual (even warns for classes without virtual functions)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd1418 -wd1419 -wd383 -wd981 -wd869")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -wd11074 -wd11076 -wd444")
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfloat-equal -Wshadow")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wswitch-default -Wconversion")
endif()


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

