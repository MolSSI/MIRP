add_library(mirp
               math.c
               gpt.c
               mpfr_help.c
               arb_help.c

               boys_double.c
               boys_mp.c
               boys_interval.c

               eri.c
               eri_mp.c
)

target_link_libraries(mirp flint arb mpfr gmp m)
