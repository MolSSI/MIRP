
macro(mirp_find_library libname)
    find_library(MIRP_${libname}_LIBRARY ${libname})

    if(NOT MIRP_${libname}_LIBRARY)
        message(FATAL_ERROR "Library ${libname} not found!")
    endif()

    get_filename_component(MIRP_${libname}_DIR ${MIRP_${libname}_LIBRARY} DIRECTORY)
    get_filename_component(MIRP_${libname}_DIR ${MIRP_${libname}_DIR} DIRECTORY)
    message(STATUS "Found ${libname} library: ${MIRP_${libname}_LIBRARY}")
    message(STATUS "Found ${libname} path: ${MIRP_${libname}_DIR}")
    list(APPEND MIRP_DEPS_LIBRARIES "${MIRP_${libname}_LIBRARY}")
    list(APPEND MIRP_DEPS_INCLUDE_DIRS "${MIRP_${libname}_DIR}/include")
endmacro()

mirp_find_library(arb)
mirp_find_library(flint)
mirp_find_library(mpfr)
mirp_find_library(gmp)
