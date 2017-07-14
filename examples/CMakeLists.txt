add_executable(simple_eri simple_eri.c)
target_link_libraries(simple_eri mirp
                                 ${PYTHON_LIBRARIES}
                                 ${MIRP_arb_LIBRARY}
                                 ${MIRP_flint_LIBRARY}
                                 ${MIRP_mpfr_LIBRARY}
                                 ${MIRP_gmp_LIBRARY}
)

