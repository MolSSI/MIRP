/*! \file
 *
 * \brief Functions related to gaussians and shells
 */

#include "mirp/shell.h"
#include "mirp/math.h"
#include <math.h>


int mirp_iterate_gaussian(int * lmn)
{
    const int am = lmn[0] + lmn[1] + lmn[2];
    if(lmn[2] >= am)
        return 0;

    if(lmn[2] < (am - lmn[0]))
    {
        lmn[1]--;
        lmn[2]++;
    }
    else
    {
        lmn[0]--;
        lmn[1] = am-lmn[0];
        lmn[2] = 0;
    }
    return 1;
}


void mirp_gaussian_fill_lmn(int am, int * lmn)
{
    const long ncart = MIRP_NCART(am);

    int g[3] = {am, 0, 0};
    for(long i = 0; i < ncart; i++)
    {
        lmn[i*3+0] = g[0];
        lmn[i*3+1] = g[1];
        lmn[i*3+2] = g[2];
        mirp_iterate_gaussian(g);
    }
}


void mirp_normalize_shell(int am, int nprim, int ngeneral,
                          arb_srcptr alpha,
                          arb_srcptr coeff,
                          arb_ptr coeff_out,
                          slong working_prec)
{
    arb_t tmp1, tmp2;
    arb_t norm_fac, sum;
    arb_t m, m2;

    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(norm_fac);
    arb_init(sum);
    arb_init(m);
    arb_init(m2);

    /* m = am + 1.5 */
    arb_set_si(m, am);
    arb_set_si(tmp1, 3);
    arb_div_si(tmp1, tmp1, 2, working_prec);
    arb_add(m, m, tmp1, working_prec);

    /* m2 = 0.5 * m */
    arb_div_si(m2, m, 2, working_prec);

    /* Normalization factor
     *
     * norm_fac = pi^(3/2) * (2l-1)!! / 2^l
     */
    arb_const_sqrt_pi(norm_fac, working_prec);
    arb_pow_ui(norm_fac, norm_fac, 3, working_prec);
    for(int i = 1; i <= am; i++)
    {
        arb_set_si(tmp1, 2*i-1);
        arb_div_si(tmp1, tmp1, 2, working_prec);
        arb_mul(norm_fac, norm_fac, tmp1, working_prec);
    }

    for(int n = 0; n < ngeneral; n++)
    {
        arb_zero(sum);

        for(int i = 0; i < nprim; i++)
        {
            for(int j = 0; j < nprim; j++)
            {
                /* tmp1 = coeff1 * coeff2 */
                arb_mul(tmp1, coeff+(n*nprim+i), coeff+(n*nprim+j), working_prec);
                /* tmp2 = pow(alpha1 * alpha2, m2) */
                arb_mul(tmp2, alpha+i, alpha+j, working_prec);
                arb_pow(tmp2, tmp2, m2, working_prec);
                arb_mul(tmp1, tmp1, tmp2, working_prec);

                /* tmp2 = pow(alpha2+alpha2, m) */
                arb_add(tmp2, alpha+i, alpha+j, working_prec);
                arb_pow(tmp2, tmp2, m, working_prec);
                arb_div(tmp1, tmp1, tmp2, working_prec);

                /* accumulate in sum */
                arb_add(sum, sum, tmp1, working_prec);
            }
        }

        /* tmp1 = 1.0 / sqrt(sum * norm_fac) */
        arb_mul(tmp1, sum, norm_fac, working_prec);
        arb_sqrt(tmp1, tmp1, working_prec);
        arb_set_si(tmp2, 1);
        arb_div(tmp1, tmp2, tmp1, working_prec);

        for (int i = 0; i < nprim; ++i)
        {
            arb_pow(tmp2, alpha+i, m2, working_prec);
            arb_mul(tmp2, tmp2, tmp1, working_prec);
            arb_mul(coeff_out+(n*nprim+i), coeff+(n*nprim+i), tmp2, working_prec);
        }
    }


    arb_clear(tmp1);
    arb_clear(tmp2);
    arb_clear(norm_fac);
    arb_clear(sum);
    arb_clear(m);
    arb_clear(m2);
}

#ifdef __cplusplus
}
#endif

