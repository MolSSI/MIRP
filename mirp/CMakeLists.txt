
list(APPEND MIRP_FILELIST
               math.c
               gpt.c
               shell.c
               loops_double.c
               loops_interval.c

               kernels/integral4_wrappers.c

               kernels/boys_double.c
               kernels/boys_interval.c

               kernels/eri_double.c
               kernels/eri_interval.c
)

add_library(mirp SHARED ${MIRP_FILELIST})
target_link_libraries(mirp PRIVATE
                           ${MIRP_arb_LIBRARY}
                           ${MIRP_flint_LIBRARY}
                           ${MIRP_mpfr_LIBRARY}
                           ${MIRP_gmp_LIBRARY}
)
