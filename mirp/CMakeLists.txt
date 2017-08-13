
list(APPEND MIRP_FILELIST
               math.c
               gpt.c
               shell.c
               loops_d.c
               loops.c

               kernels/integral4_wrappers.c

               kernels/boys_d.c
               kernels/boys.c

               kernels/eri_d.c
               kernels/eri.c
)

add_library(mirp SHARED ${MIRP_FILELIST})
target_link_libraries(mirp PRIVATE
                           ${MIRP_arb_LIBRARY}
                           ${MIRP_flint_LIBRARY}
                           ${MIRP_mpfr_LIBRARY}
                           ${MIRP_gmp_LIBRARY}
)
