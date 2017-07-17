add_library(test_common OBJECT cmdline.cpp
                               file_io.cpp
                               boys_test.cpp
)

add_executable(mirp_run_test mirp_run_test.cpp
                             $<TARGET_OBJECTS:test_common>
)

add_executable(mirp_create_test mirp_create_test.cpp
                                $<TARGET_OBJECTS:test_common>
)

add_executable(mirp_play mirp_play.cpp
                         $<TARGET_OBJECTS:test_common>
)

target_link_libraries(mirp_run_test mirp
                                ${PYTHON_LIBRARIES}
                                ${MIRP_arb_LIBRARY}
                                ${MIRP_flint_LIBRARY}
                                ${MIRP_mpfr_LIBRARY}
                                ${MIRP_gmp_LIBRARY}
)

target_link_libraries(mirp_create_test mirp
                                       ${PYTHON_LIBRARIES}
                                       ${MIRP_arb_LIBRARY}
                                       ${MIRP_flint_LIBRARY}
                                       ${MIRP_mpfr_LIBRARY}
                                       ${MIRP_gmp_LIBRARY}
)

target_link_libraries(mirp_play mirp
                                ${PYTHON_LIBRARIES}
                                ${MIRP_arb_LIBRARY}
                                ${MIRP_flint_LIBRARY}
                                ${MIRP_mpfr_LIBRARY}
                                ${MIRP_gmp_LIBRARY}
)
