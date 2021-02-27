
macro(mirp_find_library libname header)
    find_library(MIRP_${libname}_LIBRARY ${libname})

    if(NOT MIRP_${libname}_LIBRARY)
        message(FATAL_ERROR "Library ${libname} not found!")
    endif()

    get_filename_component(MIRP_${libname}_DIR ${MIRP_${libname}_LIBRARY} DIRECTORY)
    get_filename_component(MIRP_${libname}_GUESSED_INCPATH ${MIRP_${libname}_DIR}/../include DIRECTORY)

    find_path(MIRP_${libname}_INCLUDE_PATH
        NAMES ${header}
        HINTS ${MIRP_${libname}_GUESSED_INCPATH}
        DOC "MIRP include path for ${libname}")

    if(NOT MIRP_${libname}_INCLUDE_PATH)
        message(FATAL_ERROR "Include ${header} for ${libname} not found!  Please set MIRP_${libname}_INCLUDE_PATH to the correct include path.")
    endif()

    message(STATUS "Found ${libname} library: ${MIRP_${libname}_LIBRARY}")
    message(STATUS "Found ${libname} include path: ${MIRP_${libname}_INCLUDE_PATH}")

    add_library(mirp_imported_${libname} SHARED IMPORTED GLOBAL)
    set_target_properties(mirp_imported_${libname} PROPERTIES IMPORTED_LOCATION ${MIRP_${libname}_LIBRARY})
    set_target_properties(mirp_imported_${libname} PROPERTIES INTERFACE_LINK_LIBRARIES ${MIRP_${libname}_LIBRARY})
    set_target_properties(mirp_imported_${libname} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${MIRP_${libname}_INCLUDE_PATH})
    set_target_properties(mirp_imported_${libname} PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${MIRP_${libname}_INCLUDE_PATH})
    list(APPEND MIRP_DEPS_TARGETS "mirp_imported_${libname}")

endmacro()

mirp_find_library(arb arb.h)
mirp_find_library(flint flint/flint.h)
mirp_find_library(mpfr mpfr.h)
mirp_find_library(gmp gmp.h)


