#include "mirp/boys.h"
#include "mirp/arb_help.h"
#include "mirp/mpfr_help.h"
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    const int prec = 256;
    const int m = 10;

    double t, F[m+1];
    arb_t  t_arb, F_arb[m+1];
    mpfr_t t_mp, F_mp[m+1];

    arb_init(t_arb);
    mirp_init_arb_arr(F_arb, m+1);
    
    mpfr_init2(t_mp, prec);
    mirp_init_mpfr_arr(F_mp, m+1, prec);

    for(unsigned long i = 0; i < 1e4; i++)
    {
        t = (double)i / 10.0;

        arb_set_ui(t_arb, i);
        arb_div_ui(t_arb, t_arb, 10, prec);

        mpfr_set_ui(t_mp, i, MPFR_RNDN);
        mpfr_div_ui(t_mp, t_mp, 10, MPFR_RNDN);

        mirp_boys_double(F, m, t);
        mirp_boys_mp(F_mp, m, t_mp, prec);
        mirp_boys_interval(F_arb, m, t_arb, prec);

        for(int j = 0; j <= m; j++)
        {
            char * astr = arb_get_str(F_arb[j], 128, 0);
            printf("-------------------------------------\n");
            printf("i = %ld\n", i);
            printf("n = %d\n", j);
            printf(     "           double precision:  %-26.16e\n", F[j]); 
            mpfr_printf("    multi-precision [%5d]:  %Re\n", prec, F_mp[j]);
            printf(     "interval arithmetic [%5d]: %s\n", prec, astr);
            printf("        bits of accuracy: %ld\n", arb_rel_accuracy_bits(F_arb[j]));
            printf("-------------------------------------\n");
            free(astr);
        }
    }

    arb_clear(t_arb);    
    mpfr_clear(t_mp);
    mirp_clear_arb_arr(F_arb, m+1);
    mirp_clear_mpfr_arr(F_mp, m+1);
    
    return 0;
}
