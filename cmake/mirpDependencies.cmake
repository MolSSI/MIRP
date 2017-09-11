
macro(mirp_find_library libname)
    find_library(MIRP_${libname}_LIBRARY ${libname})

    if(NOT MIRP_${libname}_LIBRARY)
        message(FATAL_ERROR "Library ${libname} not found!")
    endif()

    get_filename_component(MIRP_${libname}_DIR ${MIRP_${libname}_LIBRARY} DIRECTORY)
    get_filename_component(MIRP_${libname}_DIR ${MIRP_${libname}_DIR} DIRECTORY)
    message(STATUS "Found ${libname} library: ${MIRP_${libname}_LIBRARY}")
    message(STATUS "Found ${libname} path: ${MIRP_${libname}_DIR}")

    add_library(mirp_imported_${libname} SHARED IMPORTED GLOBAL)
    set_target_properties(mirp_imported_${libname} PROPERTIES IMPORTED_LOCATION ${MIRP_${libname}_LIBRARY})
    set_target_properties(mirp_imported_${libname} PROPERTIES INTERFACE_LINK_LIBRARIES ${MIRP_${libname}_LIBRARY})
    set_target_properties(mirp_imported_${libname} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${MIRP_${libname}_DIR}/include)
    set_target_properties(mirp_imported_${libname} PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${MIRP_${libname}_DIR}/include)
    list(APPEND MIRP_DEPS_TARGETS "mirp_imported_${libname}")

endmacro()

mirp_find_library(arb)
mirp_find_library(flint)
mirp_find_library(mpfr)
mirp_find_library(gmp)


