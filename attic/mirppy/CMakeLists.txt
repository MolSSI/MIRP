find_package(PythonInterp 3.5 REQUIRED)
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(mirppy MODULE 
                   export_boys.c
                   export_eri.c
                   export.c)

set_target_properties(mirppy PROPERTIES PREFIX "")
target_link_libraries(mirppy mirp ${PYTHON_LIBRARIES}
                                  ${MIRP_arb_LIBRARY}
                                  ${MIRP_flint_LIBRARY}
                                  ${MIRP_mpfr_LIBRARY}
                                  ${MIRP_gmp_LIBRARY}
)
