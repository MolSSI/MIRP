
list(APPEND MIRP_FILELIST
               math.c
               gpt.c
               mpfr_help.c
               arb_help.c
               shell.c

               kernels/boys_double.c
               kernels/boys_mp.c
               kernels/boys_interval.c

               kernels/eri_double.c
               kernels/eri_mp.c
               kernels/eri_interval.c
)

add_library(mirp ${MIRP_FILELIST})
