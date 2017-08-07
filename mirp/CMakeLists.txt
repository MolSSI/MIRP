
list(APPEND MIRP_FILELIST
               math.c
               gpt.c
               arb_help.c
               shell.c
               loops_double.c
               loops_interval.c

               kernels/boys_double.c
               kernels/boys_interval.c

               kernels/eri_double.c
               kernels/eri_interval.c
)

add_library(mirp SHARED ${MIRP_FILELIST})
