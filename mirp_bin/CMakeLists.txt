add_library(test_common OBJECT test_files.cpp cmdline.cpp)

add_executable(mirp_test mirp_test.cpp
                         run_boys_test.cpp
                         $<TARGET_OBJECTS:test_common>
)

target_link_libraries(mirp_test mirp
                                ${PYTHON_LIBRARIES}
                                ${MIRP_arb_LIBRARY}
                                ${MIRP_flint_LIBRARY}
                                ${MIRP_mpfr_LIBRARY}
                                ${MIRP_gmp_LIBRARY}
)
